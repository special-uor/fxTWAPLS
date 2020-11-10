#' Create hexagonal logo 
#' 
#' Create hexagonal logo for the package.
#'
#' @param subplot Image to use as the main logo.
#' @param dpi Plot resolution (dots-per-inch).
#' @param h_color Colour for hexagon border.
#' @param h_fill Colour to fill hexagon.
#' @param output Output file (hexagonal logo).
#' @param package Title for logo (package name).
#' @param p_color Colour for package name.
#' @param url URL for package repository or website.
#' @param u_size Text size for URL.
#'
#' @return Hexagonal logo.
#' @keywords internal
hex_logo <- function(subplot = system.file("images/cave-painting.png", 
                                           package = "fxTWAPLS"),
                     dpi = 600,
                     h_color = "#000000",
                     h_fill = "#696969",
                     output = system.file("images/logo.png", 
                                          package = "fxTWAPLS"),
                     package = "fxTWAPLS",
                     p_color = "#eeeeee",
                     url = "https://github.com/special-uor/fxTWAPLS",
                     u_size = 1.25) {
  hexSticker::sticker(subplot = subplot, package = package,
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = url,
                      u_angle = 30, u_color = p_color, u_size = u_size,
                      filename = output)
}