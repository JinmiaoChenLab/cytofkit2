
#' Title
#'
#' @param onSever Whether the app is running on a server.
#' @param port Manually set a port, default is randomly allocated.
#'
#' @return
#' @export
#'
#' @examples cytofkit_shiny_dashboard()
cytofkit_shiny_dashboard = function(onServer = F, port = NULL){
  if(isTRUE(onServer)){
    host <- "0.0.0.0"
  }else{
    host <- "127.0.0.1"
  }
  if(!is.null(port)){
    shiny::runApp(system.file('shiny', package = 'cytofkit2'), host = host,
                  port = port)
  } else {
    shiny::runApp(system.file('shiny', package = 'cytofkit2'), host = host)
  }
}