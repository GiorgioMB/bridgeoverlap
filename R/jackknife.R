# Delete-a-group jackknife -----------------------------------------
#
# Both samples are partitioned into G contiguous blocks aligned with the
# sampling order; posterior blocks never straddle chain boundaries (G is a
# multiple of the number of chains M and each chain holds G/M blocks). Deleting
# block b from both samples leaves s1 and s2 unchanged, so the leave-block-out
# estimate is the root of W - W_b with the same responsibility function.

.split_sizes <- function(n, k) diff(round(seq(0, n, length.out = k + 1L)))

.block_ids_post <- function(chain_lengths, G) {
  per <- G %/% length(chain_lengths)
  unlist(lapply(seq_along(chain_lengths), function(m) {
    rep((m - 1L) * per + seq_len(per), times = .split_sizes(chain_lengths[m], per))
  }), use.names = FALSE)
}

.block_ids_prop <- function(S2, G) rep(seq_len(G), times = .split_sizes(S2, G))

# Cross-fitting: G / K blocks per fold. Posterior blocks are formed within the
# chains of each fold (chains are nested in folds and ordered fold by fold), and
# the proposal draws of fold k are split into the G / K blocks numbered like the
# posterior blocks of that fold, so that deleting a block removes posterior and
# proposal draws of one fold in proportion and leaves s1, s2 unchanged.
.block_ids_prop_folds <- function(fold_prop, G) {
  K <- max(fold_prop)
  per <- G %/% K
  out <- integer(length(fold_prop))
  for (k in seq_len(K)) {
    idx <- which(fold_prop == k)
    out[idx] <- (k - 1L) * per + rep(seq_len(per), times = .split_sizes(length(idx), per))
  }
  out
}

.jackknife <- function(a_post, a_prop, lr, chat, block_post, block_prop, G, tol, maxit) {
  vapply(seq_len(G), function(b) {
    .bridge_root(a_post[block_post != b], a_prop[block_prop != b], lr, c0 = chat, tol = tol, maxit = maxit)$c
  }, numeric(1))
}
