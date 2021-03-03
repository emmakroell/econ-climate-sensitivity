source("reduced_model/logistic_prcc_reduced.R",echo=TRUE)
## table_1
## prcc_res1, logistic_res1
source("full_model/logistic_and_prcc.R",echo=TRUE)
## res 2 (w/o damages & policy)
## prcc_res2, logistic_res2
## res 3 (w/ damages and policy)
## prcc_res3, logistic_res3


library(tidyverse); theme_set(theme_bw())
get_prcc <- function(i) {
    x <- get(paste0("prcc_res",i))
    ret <- (x$PRCC
        %>% rownames_to_column("param")
        %>% as_tibble()
        %>% select(param,est=original, lwr="min. c.i.", upr="max. c.i.")
        %>% mutate(across(param,~gsub("good\\.","",.)))
    )
    return(ret)
}
get_logist <- function(i) {
    x <- get(paste0("logistic_res",i))
    ret <- (broom::tidy(x,conf.int=TRUE)
        %>% select(param=term,est=estimate,lwr=conf.low,upr=conf.high)
        %>% filter(param!="(Intercept)")
    )
    return(ret)
}

tnms <- c("without climate","without damages and policy",
          "with damages and policy")
L <- map(1:3, ~bind_rows(list(PRCC=get_prcc(.),logistic=get_logist(.)), .id="type"))
names(L) <- tnms
L2 <- (bind_rows(L,.id="model"))
## order by average PRCC
PRCC_levs <- (L2
    %>% filter(type=="PRCC")
    %>% group_by(param)
    %>% summarise(est=mean(est,na.rm=TRUE))
    %>% arrange(est)
    %>% pull(param) %>% c()
)
L2 <- mutate(L2,across(param, ~factor(.,levels=PRCC_levs))) %>%
                drop_na(param)  ## 


gg0 <- (ggplot(L2,
       aes(est, param, xmin=lwr, xmax=upr)) +
    geom_pointrange() +
    facet_grid(model~type ,scale="free") +
    geom_vline(xintercept=0,lty=2) +
    labs(x="",y="") + 
    theme(panel.spacing=grid::unit(0,"lines"),
          strip.text.y.right = element_text(angle = 0))
)

#ggsave("dotwhisker.pdf",height=4,width=8)
