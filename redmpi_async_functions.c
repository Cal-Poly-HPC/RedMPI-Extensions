/* Author: Corey Ford
 * Latest version of this file available at:
 * https://gist.github.com/9022679c3d84058e48fa037d8146cf42
*/

#ifdef MRMPI_SDC_AsyncHash
/**
 * Verifies the integrity of messages between replicas using a delayed hash.
 * May be called repeatedly on the same request.
 * @param  req    The request type struct to veryify (IN).
 * @return        MPI_SUCCESS for success (including inconclusive situations
 *                where the procedure should be called again), or MPI_ERR_OTHER
 *                for error.
 */
int MRMPI_Reqs_verifyintegrity_AsyncHash(MRMPI_ReqTypePtr req) {
  if (NULL == req) {
    MRMPI_LOG_ERROR(The req parameter is null)
    return MPI_ERR_OTHER;
  }

  int hashlen = gcry_md_get_algo_dlen(GCRY_MD_SHA1);
  int error, flag;

  if (!req->myhash && req->checksdc) {
    if (NULL == (req->myhash = malloc(hashlen))) {
      MRMPI_LOG_ERROR(Unable to allocate the hash storage)
      return MPI_ERR_OTHER;
    }

    gcry_md_hash_buffer(GCRY_MD_SHA1, req->myhash, req->userbuffer,
                        req->bufsize);
  }

  if (MPI_SUCCESS != (error = PMPI_Test(&(req->reqsarray[1]), &flag,
                                        MPI_STATUS_IGNORE))) {
    MRMPI_LOG_MPI_ERROR(error,Unable to test async hash request)
    return error;
  }

  if (!flag) {
    req->incomplete = 1;
    return MPI_SUCCESS;
  }

  if (req->checksdc) {
    if (memcmp(req->hash, req->myhash, hashlen)) {
      MRMPI_LOG_ERROR(Hashes do not match)
      return MPI_ERR_OTHER;
    }
  }

  MRMPI_Reqs_free(&req);
  return MPI_SUCCESS;
}

int MRMPI_Reqs_check_AsyncHash() {
  int error, i;

  for (i = 0; i < mrmpi_reqs.count; i++) {
    if (mrmpi_reqs.array[i].inuse && mrmpi_reqs.array[i].incomplete) {
      /* Check for consistency, possibly initiating correction. */
      if (MPI_SUCCESS != (error = MRMPI_Reqs_verifyintegrity_AsyncHash(&mrmpi_reqs.array[i]))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to verify integrity)
        return error;
      }
    }
  }

  return MPI_SUCCESS;
}

int MRMPI_Reqs_finalize_AsyncHash() {
  int error, i;

  for (i = 0; i < mrmpi_reqs.count; i++) {
    if (mrmpi_reqs.array[i].inuse && mrmpi_reqs.array[i].incomplete) {
      /* Wait for all requests. */
      if (MPI_SUCCESS != (error = PMPI_Waitall(MRMPI_NUM_REQS, mrmpi_reqs.array[i].reqsarray,
                                               MPI_STATUSES_IGNORE))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to wait for all requests to complete)
        return error;
      }
      /* Check for consistency, possibly initiating correction. */
      if (MPI_SUCCESS != (error = MRMPI_Reqs_verifyintegrity_AsyncHash(&mrmpi_reqs.array[i]))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to verify integrity)
        return error;
      }
    }
  }

  return MPI_SUCCESS;
}
#endif

#ifdef MRMPI_SDC_Async
/**
 * Verifies the integrity of messages between replicas. May be called
 * repeatedly on the same request.
 * @param  req    The request type struct to veryify (IN).
 * @return        MPI_SUCCESS for success (including inconclusive situations
 *                where the procedure should be called again), or MPI_ERR_OTHER
 *                for error.
 */
int MRMPI_Reqs_verifyintegrity_Async (MRMPI_ReqTypePtr req) {
  if (NULL == req) {
    MRMPI_LOG_ERROR(The req parameter is null)
    return MPI_ERR_OTHER;
  }
  if (0 == (*req).bufsize || NULL == (*req).buffers) {
    /* Return success if the buffer size is 0 which implies no message data. */
    /* Return success if there are no buffers to verify. */

#ifdef MRMPI_SDC_Async_CopySend
    /* Clean up a send request only if all underlying sends are complete. */
    if (req->sendbuffer) {
      int error, flag;

      if (MPI_SUCCESS != (error = PMPI_Testall(mrmpi.idegree, req->reqsarray,
                                               &flag, MPI_STATUSES_IGNORE))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to testall async send)
        return error;
      }

      if (!flag) {
        req->incomplete = 1;
        return MPI_SUCCESS;
      }
    }
#endif
    MRMPI_Reqs_free(&req);
    return MPI_SUCCESS;
  }

  int hashlen = gcry_md_get_algo_dlen(GCRY_MD_SHA1);
  int total = 0, matching = 0, i;

  /* Check status of all requests. If the user's buffer has not yet been filled
   * and some request has completed, fill the user's buffer from that of the
   * first such request encountered. Count how many requests have a buffer
   * matching this first request.
   */
  for (i = 0; i < mrmpi.idegree; i++) {
    if (!req->hashes[i]) {
      int error, flag;

      if(MPI_SUCCESS != (error = PMPI_Test(&(req->reqsarray[i]), &flag,
                                           MPI_STATUS_IGNORE))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to test async request)
        return error;
      }

      if (flag) {
        if (NULL == (req->hashes[i] = malloc(hashlen))) {
          MRMPI_LOG_ERROR(Unable to allocate the hash storage)
          return MPI_ERR_OTHER;
        }
        gcry_md_hash_buffer(GCRY_MD_SHA1, req->hashes[i], req->buffers[i],
                            req->bufsize);
      }
    }

    if (req->hashes[i]) {
      total++;
      if (req->firstindex < 0) {
        memcpy(req->userbuffer, req->buffers[i], req->bufsize);
        req->firstindex = i;
        matching = 1;
      } else if (0 == memcmp(req->hashes[i], req->hashes[req->firstindex],
                             hashlen)) {
        matching++;
      }

      /* Buffer contents are no longer needed. */
      if (req->buffers[i]) {
        free(req->buffers[i]);
        req->buffers[i] = NULL;
      }
    }
  }

  /* If all messages have been received, discard the request entry.
   * The buffers cannot be freed before this time.
   */
  if (total == mrmpi.idegree) {
    MRMPI_Reqs_free(&req);

    /* If a majority of messages match the first, success.
     * Otherwise, initiate correction.
     */
    if (matching > mrmpi.idegree / 2) {
      return MPI_SUCCESS;
    }
    return MPI_ERR_OTHER;
  }

  /* Insufficient information, mark request as incomplete. */
  req->incomplete = 1;
  return MPI_SUCCESS;
}

int MRMPI_Reqs_check_Async() {
  int error, i;

  for (i = 0; i < mrmpi_reqs.count; i++) {
    if (mrmpi_reqs.array[i].inuse && mrmpi_reqs.array[i].incomplete) {
      /* Check for consistency, possibly initiating correction. */
      if (MPI_SUCCESS != (error = MRMPI_Reqs_verifyintegrity_Async(&mrmpi_reqs.array[i]))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to verify integrity)
        return error;
      }
    }
  }

  return MPI_SUCCESS;
}

int MRMPI_Reqs_finalize_Async() {
  int error, i;

  for (i = 0; i < mrmpi_reqs.count; i++) {
    if (mrmpi_reqs.array[i].inuse && mrmpi_reqs.array[i].incomplete) {
      /* Wait for all requests. */
      if (MPI_SUCCESS != (error = PMPI_Waitall(MRMPI_NUM_REQS, mrmpi_reqs.array[i].reqsarray,
                                               MPI_STATUSES_IGNORE))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to wait for all requests to complete)
        return error;
      }
      /* Check for consistency, possibly initiating correction. */
      if (MPI_SUCCESS != (error = MRMPI_Reqs_verifyintegrity_Async(&mrmpi_reqs.array[i]))) {
        MRMPI_LOG_MPI_ERROR(error,Unable to verify integrity)
        return error;
      }
    }
  }

  return MPI_SUCCESS;
}
#endif

