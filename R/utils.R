#' Create hexagonal logo for the package
#'
#' @param subplot image to use as the main logo
#' @param dpi plot resolution (dots-per-inch)
#' @param h_color color for hexagon border
#' @param h_fill color to fill hexagon
#' @param output output file (hexagonal logo)
#' @param package title for logo (package name)
#' @param p_color color for package name
#' @param url URL for package repository or website
#' @param u_size text size for URL
#'
#' @return hexagonal logo
#' @export
#'
#' @examples
#' hex_logo()
#' \dontrun{
#' hex_logo("inst/images/cave-painting.png", output = "inst/images/logo.png")
#' }
hex_logo <- function(subplot = system.file("images/cave-painting.png", 
                                           package = "reconbio"),
                     dpi = 600,
                     h_color = "#000000",
                     h_fill = "#696969",
                     output = system.file("images/logo.png", package = "reconbio"),
                     package = "reconbio",
                     p_color = "#eeeeee",
                     url = "https://github.com/special-uor/reconbio",
                     u_size = 1.30) {
  hexSticker::sticker(subplot = subplot, package = package,
                      h_color = h_color,  h_fill = h_fill,
                      dpi = dpi,
                      s_x = 1.0, s_y = .85, s_width = .5,
                      p_x = 1.0, p_y = 1.52, p_size = 6, p_color = p_color,
                      url = url,
                      u_angle = 30, u_color = p_color, u_size = u_size,
                      filename = output)
}
