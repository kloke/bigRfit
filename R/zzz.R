setGeneric(name="getScoresNS", def=function(object, x) {standardGeneric("getScoresNS")})
setGeneric(name="getScoresDerivNS", def=function(object, x) {standardGeneric("getScoresDerivNS")})

setMethod("getScoresNS", "scores", 
  function(object, x) {
    if(is.null(object@param)) {
      a<-object@phi(x)
     } else {
      a<-object@phi(x, object@param)
     }
     a
   }
)


setMethod("getScoresDerivNS", "scores",
  function(object, x) {
    if(is.null(object@param)) {
      aP<-object@Dphi(x)
     } else {
      aP<-object@Dphi(x, object@param)
     }
     aP
   }
)

