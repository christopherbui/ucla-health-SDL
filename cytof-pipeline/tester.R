library(ggplot2)

df <- data.frame(x = 1:10, y = rnorm(10))

ggplot(df, aes(x, y)) + geom_point()
