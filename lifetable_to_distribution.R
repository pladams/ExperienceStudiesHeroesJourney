lifetable_to_distribution <- function(qx_vector,
                                      terminate_table=TRUE) {
  event_probs <- data.table(
    time=1:length(qx_vector),
    qx=qx_vector
  )
  event_probs[,
              tpx:=cumprod(1-qx)
  ]
  event_probs[,
              q_xpt:=qx*shift(tpx,fill=1)
  ]
  
  
  if(terminate_table) {
    rbind(
      event_probs,
      data.table(
        time=length(qx_vector)+1,
        qx = 1,
        tpx=0,
        q_xpt = 1-event_probs[,sum(q_xpt)]
      )
    ) -> event_probs
  }
  
  return(event_probs)
}
