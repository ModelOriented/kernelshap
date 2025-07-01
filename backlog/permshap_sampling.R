library(ranger)
set.seed(1)
object <- lm(Sepal.Length ~ ., data = iris)
pred_fun <- predict

object <- ranger(Sepal.Length ~ ., data = iris)


system.time(permshap(object, X = iris[-1], bg_X = iris[, -1]))
system.time(permshap(object, X = iris[-1], bg_X = iris[, -1], exact = F))

# Sepal.Width Petal.Length Petal.Width   Species
#   0.2195135    -1.955357   0.3149451 0.5823533

# Sepal.Width Petal.Length Petal.Width    Species
#   0.1734389   -0.5202173  -0.2416169 -0.1420209
