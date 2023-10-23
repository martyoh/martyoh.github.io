simple_plot <- function(start, stop) {
  df <- data.frame(x=start:stop, y=start:stop, z=stop:start)
  print(plot(df$x, df$y))
}

