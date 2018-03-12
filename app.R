
library(shiny)

library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

plot_res <- function(res, b) {

  res$pi <- tail(res$pi, nrow(res$pi) - b)
  res$mu <- tail(res$mu, nrow(res$mu) - b)
  res$sd <- tail(res$sd, nrow(res$sd) - b)
  res$ll <- tail(res$ll, length(res$ll) - b)

  p1 <- res$pi %>% as.data.frame() %>% mutate(iter = row_number()) %>% gather(sample, mix, -iter) %>%
    ggplot() +
    geom_line(aes(iter, mix, color = sample)) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = 'none')
  p2 <- res$mu %>% as.data.frame() %>% mutate(iter = row_number()) %>% gather(sample, mean, -iter) %>%
    ggplot() +
    geom_line(aes(iter, mean, color = sample)) +
    scale_color_brewer(palette = "Dark2") +
    scale_y_continuous(limits = c(0,1)) +
    theme(legend.position = 'none')
  p3 <- res$sd %>% as.data.frame() %>% mutate(iteration = row_number()) %>% gather(sample, StdDev, -iteration) %>%
    ggplot() +
    geom_line(aes(iteration, StdDev, color = sample)) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = 'none')

  p4 <- res$pi %>% as.data.frame() %>% gather(sample, val) %>%
    ggplot() +
    geom_density(aes(val, color = sample)) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(0,1)) +
    coord_flip()

  p5 <- res$mu %>% as.data.frame() %>% gather(sample, val) %>%
    ggplot() +
    geom_density(aes(val, color = sample)) +
    scale_color_brewer(palette = "Dark2") +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    scale_x_continuous(limits = c(0,1)) +
    coord_flip()

  p6 <- res$ll %>% as.data.frame() %>% mutate(iteration = row_number()) %>% gather(sample, logLik, -iteration) %>%
    ggplot() +
    geom_line(aes(iteration, logLik)) +
    theme(legend.position = 'none')

  p7 <- res$ll %>% as.data.frame() %>% gather(sample, val) %>%
    ggplot() +
    geom_density(aes(val)) +
    theme(legend.position = 'none',
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank()) +
    coord_flip()

  g1 <- plot_grid(p6 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                  p1 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                  p2 + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
                  p3, nrow = 4, align = 'v')

  g2 <- plot_grid(p7, p4, p5, ggplot(), nrow = 4, align = 'v')

  final <- plot_grid(g1, g2, nrow = 1, rel_widths = c(8,2))

  return(final)

}

build_tbl <- function(res, b) {
  res$pi <- tail(res$pi, nrow(res$pi) - b)
  res$mu <- tail(res$mu, nrow(res$mu) - b)
  res$sd <- tail(res$sd, nrow(res$sd) - b)
  res$ll <- tail(res$ll, length(res$ll) - b)

  myTbl <- data.frame(
    "Percentile" = c(2.5, 50, 97.5),
    "logLikelihood" = quantile(res$ll, probs = c(0.025,0.5,0.975), na.rm = T),
    row.names = c("2.50% Perc.", "50.0% Perc.", "97.5% Perc.")
  )
  for (i in c(1:ncol(res$pi))) {
    myTbl[[paste("mix", i)]]  <- quantile(res$pi[,i], probs = c(0.025,0.5,0.975), na.rm = T)
    myTbl[[paste("mean", i)]] <- quantile(res$mu[,i], probs = c(0.025,0.5,0.975), na.rm = T)
  }
  return(myTbl)
}

gibbs <- function(x,k,niter=1000,muprior = list(mean=0.5,var=1),sdprior = list(n=2,var=0.01)){
  pi = rep(1/k,k)
  mu = rnorm(k,muprior$mean,sqrt(muprior$var))
  sd = sqrt(1/rgamma(k,sdprior$n/2,sdprior$var*sdprior$n/2))
  z = sample_z(x,pi,mu,sd)
  res = list(mu = matrix(nrow=niter, ncol=k), sd = matrix(nrow=niter, ncol=k),
             pi = matrix(nrow=niter,ncol=k), z = matrix(nrow=niter, ncol=length(x)),
             ll = rep(NA, niter))
  res$mu[1,]=mu
  res$sd[1,]=sd
  res$pi[1,]=pi
  res$z[1,]=z
  for(i in 2:niter){
    pi = sample_pi(z,k)
    mu = sample_mu(x,z,k,sd,muprior)
    sd = sample_sd(x,z,k,mu,sdprior)
    z = sample_z(x,pi,mu,sd)
    res$mu[i,] = mu
    res$sd[i,] = sd
    res$pi[i,] = pi
    res$z[i,] = z
    res$ll[i] = calc_loglik(x, pi, mu, sd)
  }
  return(res)
}

normalize <- function(x){return(x/sum(x))}

sample_z <- function(x,pi,mu,sd){
  p.x.given.z = matrix(ncol = length(pi), nrow = length(x))
  for(i in 1:length(pi)) {
    p.x.given.z[,i] = pi[i] * dnorm(x, mu[i], sd[i])
  }
  p.x.given.z = apply(p.x.given.z,1,normalize)
  if (!is.vector(p.x.given.z)) {
    p.x.given.z = t(p.x.given.z)
  }
  z = rep(0, length(x))
  for(i in 1:length(z)){
    z[i] = sample(1:length(pi), size=1,prob=as.matrix(p.x.given.z)[i,],replace=TRUE)
  }
  return(z)
}

sample_pi <- function(z,k){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi = gtools::rdirichlet(1,counts+1)
  return(pi)
}

sample_mu <- function(x, z, k, sd, prior){
  df = data.frame(x=x,z=z)
  mu = rep(0,k)
  for(i in 1:k){
    sample.size = sum(z==i)
    sample.mean = ifelse(sample.size==0,0,mean(x[z==i]))

    post.mean = (prior$mean/prior$var + (sample.size*sample.mean)/sd[i]^2)/(1/prior$var + sample.size/sd[i]^2)
    post.var = 1/(1/prior$var + sample.size/sd[i]^2)

    mu[i] = rnorm(1,post.mean,sqrt(post.var))
  }
  return(mu)
}

sample_sd <- function(x, z, k, mu, prior){
  df = data.frame(x=x,z=z)
  sd = rep(0,k)
  for(i in 1:k){
    sample.size = sum(z==i)
    sample.mean = ifelse(sample.size==0,0,mean(x[z==i]))
    sample.var = ifelse(sample.size<2,0,var(x[z==i]))

    post.n = prior$n + sample.size
    post.var = (prior$n*prior$var + (sample.size-1)*sample.var + sample.size * (sample.mean - mu[i])^2) / post.n

    sd[i] = sqrt(1/rgamma(1, post.n/2, post.n*post.var/2))
  }
  return(sd)
}

calc_loglik <- function(x, pi, mu, sd){
  return( sum( log (
    sapply(x, function(i, pi, mu, sd) {
      sum(pi * dnorm(i, mu, sd))
    }, pi = pi, mu = mu, sd = sd)
  )))
}

mixture.EM <- function(L, X, d, f = c(1,1,1)) {

  # store log-likehoods for each iteration
  log_liks <- c()
  ll       <- calc_loglik(X, L[1,], L[2,], L[3,])
  log_liks <- c(log_liks, ll)
  delta.ll <- 100000

  while(delta.ll > d) {
    L        <- EM.iter(X, L, f)
    ll       <- calc_loglik(X, L[1,], L[2,], L[3,])
    log_liks <- c(log_liks, ll)
    delta.ll <- log_liks[length(log_liks)]  - log_liks[length(log_liks)-1]

    if (is.na(delta.ll)) {
      delta.ll <- 100000
      L[2,] <- runif(ncol(L), 0.2, 0.8)
    }
  }
  return(list(L, log_liks))
}

EM.iter <- function(X, L, f) {

  P <- L

  # E-step
  # get estimates for latent variable posterior based on current parameters
  gamma <- apply(L, 2, function(p, D){
              return(p[1] * dnorm(D, mean = p[2], sd = p[3]))
           }, D = X)

  gamma = gamma / rowSums(gamma)

  # M-step
  # only estimate parameter group is corresponding f[x] is set
  if (f[1]) {
    # get new maximum for mixture proportions
    L[1,] <- colSums(gamma) / sum(gamma)
  }

  if (f[2]) {
    # get new maximum for means
    L[2,] <- 1/colSums(gamma) * colSums(gamma*X)
  }

  if (f[3]) {
    # get new maximum for standard deviation
    for (i in 1:ncol(L)) {
      L[3,i] <- sqrt( 1/sum(gamma[,i]) * sum(gamma[,i] * (X - L[2,i])^2 ))
    }
  }

  # remove NaN from L that might arise when fixing zeros
  L <- replace(L, is.na(L), P)
  return(L)
}

plot_EM <- function(res) {
  res[[2]] %>% as.data.frame() %>% mutate(iteration = row_number()) %>% gather(sample, logLik, -iteration) %>%
    ggplot() +
    geom_line(aes(iteration, logLik)) +
    theme(legend.position = 'none')

}

build_tblEM <- function(res) {

  myTbl <- res[[1]] %>% as.data.frame()
  rownames(myTbl) <- c("Mix", "Mean", "StdDev")
  return(myTbl)
}



# Define UI ----
ui <- fluidPage(

  sidebarLayout(
    sidebarPanel(
      h3("Input"),

      h4("Load bin file"),

      fileInput("infile", NULL, buttonLabel = "Load", placeholder = "or use simulation below"),

      conditionalPanel(
        condition = "output.fileUsed == true",
        numericInput("nf", h5("Maximum n from real Data"), value = 1000, min = 10, step = 10)
      ),

      conditionalPanel(
        condition = "output.fileUsed == false",

        h4("Simulation"),

        numericInput("n", h5("samples"), value = 1000, min = 10, step = 10),

        fluidRow(
          column(3,
                h5("Mix 1"),
                numericInput("a1", NULL, value = 0.33, min = 0, max = 1, step = 0.01),
                h5("Mean 1"),
                numericInput("m1", NULL, value = 0.25, min = 0, max = 1, step = 0.01),
                h5("Std. dev. 1"),
                numericInput("w1", NULL, value = 0.06, min = 0, max = 1, step = 0.01)
                ),
          column(3,
                h5("2"),
                numericInput("a2", NULL, value = 0.33, min = 0, max = 1, step = 0.01),
                h5("2"),
                numericInput("m2", NULL, value = 0.50, min = 0, max = 1, step = 0.01),
                h5("2"),
                numericInput("w2", NULL, value = 0.06, min = 0, max = 1, step = 0.01)
                ),
          column(3,
                h5("3"),
                numericInput("a3", NULL, value = 0.33, min = 0, max = 1, step = 0.01),
                h5("3"),
                numericInput("m3", NULL, value = 0.75, min = 0, max = 1, step = 0.01),
                h5("3"),
                numericInput("w3", NULL, value = 0.06, min = 0, max = 1, step = 0.01)
                )
        )
      ),

      sliderInput("breaks", h5("Breaks"), 1, 100, 50),

      selectInput("method", "Method:", c("Gibbs Sampler" = "gibbs", "EM Algorithm" = "em")),

      conditionalPanel(
        condition = "input.method == 'gibbs'",

        numericInput("cG", h5("Clusters"), value = 3, min = 1, max = 10, step = 1),
        numericInput("i", h5("Iterations"), value = 250, min = 10, step = 10),

        actionButton("submitG","Submit"),

        numericInput("b", h5("Burn-in"), value = 100, step = 10)
      ),

      conditionalPanel(
        condition = "input.method == 'em'",

        checkboxInput("fix","Run fixed model with 1, 2 or 3 clusters"),

        numericInput("cE", h5("Clusters"), value = 3, min = 1, max = 10, step = 1),
        numericInput("d", h5("Delta"), value = 1e-5),

        actionButton("submitE","Submit")

        ),

      checkboxInput("keep","Do you want to record logLikelihood?"),

      conditionalPanel(
        condition = "input.keep == true",
        fluidRow(
          column(7,
                 textInput("name", "Name of the record")
                 ),
          column(2,
                 style = "margin-top: 25px;",
                 actionButton("addN", "Add")
                 ),
          column(3,
                 style = "margin-top: 25px;",
                 actionButton("clearN", "Clear")
                 )
        )
      )


    ),
    mainPanel(
      plotOutput("histPlot"),

      conditionalPanel(
        condition = "input.method == 'gibbs'",

        plotOutput("resPlotG"),
        tableOutput("tblOutG")
      ),
      conditionalPanel(
        condition = "input.method == 'em'",

        tableOutput("tblOutE"),
        plotOutput("resPlotE")
      ),
      conditionalPanel(
        condition = "input.keep == true",
        tableOutput("tblRecE")
      ),

      uiOutput("modelHead"),
      plotOutput("modelRes")
    )
  )
)

# Define server logic ----
server <- function(input, output) {

  dyn.load('src/read_bin.so')
  binReadR <- function(file) {
    resl <- .Call("binReadR", file)
    return(unlist(resl))
  }

  output$fileUsed <- reactive({
    return(!is.null(input$infile))
  })
  outputOptions(output, 'fileUsed', suspendWhenHidden = FALSE)

  get_data <- reactive({

    if (is.null(input$infile)) {

      pi = c(input$a1, input$a2, input$a3)
      mu = c(input$m1, input$m2, input$m3)
      s = c(input$w1, input$w2, input$w3)
      z = sample(1:3,prob=pi,size=input$n,replace=TRUE)

      x = sapply(z, function(i, mu, s) {
        v = 0.0
        while (v < 0.2 || v > 0.8) {
          v = rnorm(1, mean = mu[i], sd = s[i])
        }
        return(v)
      }, mu = mu, s = s)

    } else {
      bhist <- binReadR(input$infile$datapath)
      x <- if(length(bhist) <= input$nf) bhist else sample(bhist, size = input$nf)
    }

    return(x)
  })

  output$histPlot <- renderPlot({
    x <- get_data()
    hist(x, breaks = input$breaks, xlab = "", main = "Input Histogram")
  })

  get_gibbs <- reactive({
    if (input$submitG == 0)
      return()

    isolate({
      x <- get_data()
      res = gibbs(x = x, k = input$cG, niter = input$i)
      return(res)
    })
  })

  output$resPlotG <- renderPlot({
    if (input$submitG == 0)
      return()

    res <- get_gibbs()
    print(plot_res(res, input$b))
  })

  output$tblOutG <- renderTable({
    if (input$submitG == 0)
      return()

    res <- get_gibbs()
    return(build_tbl(res, input$b))
  }, align = 'c')

  get_em <- reactive({
    if (input$submitE == 0)
      return()

    isolate({
      x <- get_data()

      if (input$fix == TRUE) {
        if (!input$cE %in% c(1,2,3))
          return(NULL)
        else if (input$cE == 1) {
          L <- matrix(c(1,
                        0.5,
                        0.01),
                        byrow = T, nrow = 3)
          res = mixture.EM(L, x, input$d, f = c(0,0,1))
        }
        else if (input$cE == 2) {
          L <- matrix(c(0.5, 0.5,
                        0.33, 0.67,
                        0.01, 0.01),
                        byrow = T, nrow = 3)
          res = mixture.EM(L, x, input$d, f = c(0,0,1))
        }
        else if (input$cE == 3) {
          L <- matrix(c(0.33, 0.33, 0.33,
                        0.25, 0.50, 0.75,
                        0.01, 0.01, 0.01),
                        byrow = T, nrow = 3)
          res = mixture.EM(L, x, input$d, f = c(0,0,1))
        }
      } else {
        if (input$cE == 3) {
          L <- matrix(c(0.20, 0.30, 0.70,
                        0.50, 0.70, 0.20,
                        0.01, 0.01, 0.01),
                        byrow = T, nrow = 3)
          res = mixture.EM(L, x, input$d)
        } else {
          L <- matrix(c(normalize(runif(input$cE, 0.2, 0.8)),
                        runif(input$cE, 0.2, 0.8),
                        rep(0.01, input$cE)),
                        byrow = T, nrow = 3)
          res = mixture.EM(L, x, input$d)
        }
      }
      return(res)
    })
  })

  output$resPlotE <- renderPlot({
    if (input$submitE == 0)
      return()

    res <- get_em()
    print(plot_EM(res))
  })

  output$tblOutE <- renderTable({
    if (input$submitE == 0)
      return()

    res <- get_em()
    return(build_tblEM(res))
  }, align = 'c', rownames = T, colnames = F)


  record <- data.frame("name" = c(NA), "logLik" = c(NA))

  observe({
    if (input$addN == 0)
      return()

    isolate({
      ll <- 0
      if (input$method == 'em') {
        res <- get_em()
        ll <- tail(res[[2]], 1)
      }
      else if (input$method == 'gibbs') {
        res <- get_gibbs()
        ll <- mean(tail(res$ll, length(res$ll) - input$b))
      }

      record[nrow(record)+1,] <<- c("name" = input$name, "logLik" = ll)

      output$tblRecE <- renderTable({
        record[complete.cases(record),]
      })
    })
  })

  observe({
    if (input$clearN == 0)
      return()

    isolate({

      record <<- data.frame("name" = c(NA), "logLik" = c(NA))

      output$tblRecE <- renderTable({
        record[complete.cases(record),]
      })
    })
  })

  observe({
    if (input$submitE == 0 && input$submitG == 0)
      return()

    input$b

    isolate({
      if (input$method == 'em') {
        res <- get_em()

        mtb <- res[[1]]

      }
      else if (input$method == 'gibbs') {
        res <- get_gibbs()

        res$pi <- tail(res$pi, nrow(res$pi) - input$b)
        res$mu <- tail(res$mu, nrow(res$mu) - input$b)
        res$sd <- tail(res$sd, nrow(res$sd) - input$b)

        mtb <- matrix(c(
          apply(res$pi, 2, mean),
          apply(res$mu, 2, mean),
          apply(res$sd, 2, mean)),
          byrow = T, nrow = 3
        )

      }

      output$modelHead <- renderUI({
        h4("Density Simulation based on Model Results:")
      })
      output$modelRes <- renderPlot({
        unlist(as.vector(apply(mtb, 2, function(i) {return(rnorm(10000*i[1], i[2], i[3]))}))) %>%
          as.data.frame() %>% setNames("simulation") %>%
          ggplot(aes(simulation)) +
          geom_density() +
         # geom_histogram(bins = 50) +
          xlim(0.2, 0.8)

      })
    })
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)
