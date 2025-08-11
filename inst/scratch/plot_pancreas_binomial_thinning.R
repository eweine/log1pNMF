out <- readr::read_rds("/Users/eweine/Documents/data/pancreas_binomial_thinning_res_final.rds")

tm_train_ll <- unlist(lapply(out, function(x){x$tm_train_ll}))
tm_test_ll <- unlist(lapply(out, function(x){x$tm_test_ll}))

train_ll_c_001 <- unlist(lapply(out, function(x){x[[as.character(0.001)]]$train_ll}))
test_ll_c_001  <- unlist(lapply(out, function(x){x[[as.character(0.001)]]$test_ll}))

train_ll_c_01 <- unlist(lapply(out, function(x){x[[as.character(0.01)]]$train_ll}))
test_ll_c_01  <- unlist(lapply(out, function(x){x[[as.character(0.01)]]$test_ll}))

train_ll_c_p1 <- unlist(lapply(out, function(x){x[[as.character(0.1)]]$train_ll}))
test_ll_c_p1  <- unlist(lapply(out, function(x){x[[as.character(0.1)]]$test_ll}))

train_ll_c_1 <- unlist(lapply(out, function(x){x[[as.character(1)]]$train_ll}))
test_ll_c_1  <- unlist(lapply(out, function(x){x[[as.character(1)]]$test_ll}))

train_ll_c_10 <- unlist(lapply(out, function(x){x[[as.character(10)]]$train_ll}))
test_ll_c_10  <- unlist(lapply(out, function(x){x[[as.character(10)]]$test_ll}))

train_ll_c_100 <- unlist(lapply(out, function(x){x[[as.character(100)]]$train_ll}))
test_ll_c_100  <- unlist(lapply(out, function(x){x[[as.character(100)]]$test_ll}))

train_ll_c_1000 <- unlist(lapply(out, function(x){x[[as.character(1000)]]$train_ll}))
test_ll_c_1000  <- unlist(lapply(out, function(x){x[[as.character(1000)]]$test_ll}))


res_df <- data.frame(
  cc = rep(c(0.001, 0.01, 0.1, 1, 10, 100, 1000, Inf), each = 18),
  loglik = c(
    train_ll_c_001,
    test_ll_c_001,
    train_ll_c_01,
    test_ll_c_01,
    train_ll_c_p1,
    test_ll_c_p1,
    train_ll_c_1,
    test_ll_c_1,
    train_ll_c_10,
    test_ll_c_10,
    train_ll_c_100,
    test_ll_c_100,
    train_ll_c_1000,
    test_ll_c_1000,
    tm_train_ll,
    tm_test_ll
  ),
  type = rep(
    rep(c("train", "test"), each = 9),
    8
  )
)

library(ggplot2)
ggplot(data = res_df, aes(x = factor(cc), y = -loglik)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~type) +
  xlab("cc") +
  ylab("Loss (-loglik)") +
  cowplot::theme_cowplot()
