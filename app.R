library(shiny)
library(deSolve)
library(shinymanager)

#cred <- data.frame(
#  user     = c("lcb",  "user"),
#  password = c("Lcb3621!",  "Lcb3621!"),
#  admin    = c(TRUE,     FALSE)
#)
#create_db(
#  credentials = cred,
#  sqlite_path = "auth.sqlite",           # ← 앱 폴더에 생김
#  passphrase  = Sys.getenv("SM_KEY")     # 환경변수로 암호화 키 보관
#)

# ODE 정의
tgi_model <- function(time, state, parameters, dosing_times, dose_amt) {
  with(as.list(c(state, parameters)), {
    # Ensure safe values to prevent division by near-zero
    V1 <- max(V1, 1e-6)
    V2 <- max(V2, 1e-6)
    V1_PL <- max(V1_PL, 1e-6)
    V2_PL <- max(V2_PL, 1e-6)
    epsil <- max(epsil, 1e-6)
    epsil_PL <- max(epsil_PL, 1e-6)
    IC50 <- max(IC50, 1e-6)
    KG <- max(KG, 1e-6)
    TAU <- max(TAU, 1e-6)
    
    
    PL <- A10 + A11 + A12
    W <- A13 + A14 + A15 + A16
    W_safe <- max(W, 1e-6)
    Rtumor <- (3 * W_safe / (4 * pi))^(1/3) / 10
    
    dA1 <- -CL / V1 * A1 - Q * (A1 / V1 - A2 / V2) - (2 * PADC * Rcap / RKrogh^2 + 6 * DADC / Rtumor^2) * (A1 / V1 - A8 / epsil) * W / 1e6
    dA2 <- Q * (A1 / V1 - A2 / V2)
    dA3 <- 0; dA4 <- 0
    dA5 <- -CL_PL / V1_PL * A5 - Q_PL * (A5 / V1_PL - A6 / V2_PL) +
      ((CL / V1 + kdis) * A7 * A1) / V1_PL -
      (2 * P_PL * Rcap / RKrogh^2 + 6 * D_PL / Rtumor^2) * (A5 / V1_PL - A12 / epsil_PL) * W / 1e6
    dA6 <- Q_PL * (A5 / V1_PL - A6 / V2_PL)
    dA7 <- -kdis * A7
    dA8 <- (2 * PADC * Rcap / RKrogh^2 + 6 * DADC / Rtumor^2) * (A1 / V1 - A8 / epsil) - Kon * A8 * (Ag - A9) / epsil + Koff * A9
    dA9 <- Kon * A8 * (Ag - A9) / epsil - Koff * A9 - Kint * A9
    dA10 <- (2 * P_PL * Rcap / RKrogh^2 + 6 * D_PL / Rtumor^2) * (A5 - A10 / epsil_PL) -
      Kint_PL * A10 + Kout_PL * A11 + kdis * A7 * (A8 + A9)
    dA11 <- Kint * A9 * A7 + Kint_PL * A10 - Kout_PL * A11 - Kon_PL * A11 * (DNA - A12) + Koff_PL * A12
    dA12 <- Kon_PL * A11 * (DNA - A12) - Koff_PL * A12
    
    growth <- (KG0 * (1 - W / M_MAX) * A13) / ((1 + ((KG0 / KG) * W)^PSI)^(1 / PSI))
    kill <- (KMAX * PL^GAMMA) / (IC50^GAMMA + PL^GAMMA)
    
    dA13 <- growth - kill * A13
    dA14 <- kill * A13 - A14 / TAU
    dA15 <- (A14 - A15) / TAU
    dA16 <- (A15 - A16) / TAU
    
    C_ADC_ng_mL <- (A1 / V1) * mw * 1e-3
    C_PL_ng_mL <- (A5 / V1_PL) * mw_PL * 1e-3
    list(c(dA1, dA2, dA3, dA4, dA5, dA6, dA7, dA8, dA9, dA10, dA11, dA12, dA13, dA14, dA15, dA16),
         W = W, PL=PL, ADC = C_ADC_ng_mL, PLconc = C_PL_ng_mL)
  })
}

# UI
ui <- secure_app(fluidPage(
  titlePanel("ADC TGI PK/PD Model"),
  div(
    style = "margin-bottom: 20px; color: #555; font-size: 14px;",
    em("Created by Jae-Yong Chung, SNUBH 2025")
  ),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(6,
               h4("Dosing"),
               numericInput("dose_amt", "Dose Amount (mg/kg)", 10),
               numericInput("dose_int", "Dosing Interval (days)", 4),
               numericInput("dose_n", "Number of Doses", 3),
               fluidRow(
                 column(6,
                        numericInput("mw", "ADC MW", 150000)
                 ),
                 column(6,
                        numericInput("mw_PL", "Payload MW", 691.7)
                 )
               )),
        column(6,
               br(), br(), br(),
               actionButton("simulate", "Simulate", class = "btn-primary", style = "float: right;")
        )
      ),
      hr(),
      fluidRow(
        column(6,
               h5("Tumor"),
               numericInput("KG0", "KG0 (exp) (1/day)", 0.08),
               numericInput("KG", "KG (lin) (mm3/day)", 225),
               numericInput("PSI", "PSI", 20),
               numericInput("M_MAX", "M_MAX (mm3)", 5000),
               numericInput("KMAX", "KMAX (1/day)", 17.6),
               numericInput("IC50", "IC50 (nM)", 399),
               numericInput("TAU", "TAU (day)", 1.21),
               numericInput("GAMMA", "GAMMA", 1)
        ),
        column(6,
               h5("PK / PD"),
               numericInput("CL", "CL", 0.039),
               numericInput("V1", "V1", 0.0478),
               numericInput("Q", "Q", 0.024),
               numericInput("V2", "V2", 0.0214),
               numericInput("kdis", "kdis", 0.24),
               numericInput("CL_PL", "CL_PL", 151.3),
               numericInput("V1_PL", "V1_PL", 1.38),
               numericInput("Q_PL", "Q_PL", 1304),
               numericInput("V2_PL", "V2_PL", 126)
        )
      ),
      hr(),
      fluidRow(
        column(6,
               h5("Transport"),
               numericInput("PADC", "PADC", 334),
               numericInput("DADC", "DADC", 0.022),
               numericInput("P_PL", "P_PL", 18144),
               numericInput("D_PL", "D_PL", 0.125),
               numericInput("Rcap", "Rcap", 8),
               numericInput("RKrogh", "RKrogh", 75)
        ),
        column(6,
               h5("Binding"),
               numericInput("epsil", "epsil", 0.24),
               numericInput("epsil_PL", "epsil_PL", 0.44),
               numericInput("Kon", "Kon", 43.2),
               numericInput("Koff", "Koff", 6.48),
               numericInput("Kint", "Kint", 199.6),
               numericInput("Kint_PL", "Kint_PL", 9.66),
               numericInput("Kon_PL", "Kon_PL", 1),
               numericInput("Koff_PL", "Koff_PL", 1),
               numericInput("Kout_PL", "Kout_PL", 1),
               numericInput("Ag", "Ag", 3.4),
               numericInput("DNA", "DNA", 196)
        )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("ADC in tumor bound", plotOutput("plotA9")),
        tabPanel("ADC in tumor free", plotOutput("plotA8")),
        tabPanel("PL in tumor interstitial", plotOutput("plotA10")),
        tabPanel("PL in cell free", plotOutput("plotA11")),
        tabPanel("PL in cell bound", plotOutput("plotA12")),
        tabPanel("PL Total", plotOutput("plotPLTotal")),
        tabPanel("Tumor Volume", plotOutput("plotW")),
        tabPanel("ADC in plasma", plotOutput("plotADC")),
        tabPanel("PL in plasma", plotOutput("plotPL")),
        tabPanel("Combined", plotOutput("plotAll"))
      )
    )
  )
)
)


# Server
server <- function(input, output,session) {
  res = secure_server(check_credentials= check_credentials("auth.sqlite",
                                                           passphrase = Sys.getenv("SM_KEY")))
  observeEvent(input$simulate, {
    parameters <- reactiveValuesToList(input)
    times <- seq(0, 50, by = 0.5)
    dosing_times <- seq(0, by = input$dose_int, length.out = input$dose_n)
    
    state <- c(A1=0, A2=0, A3=0, A4=0, A5=0, A6=0,
               A7=4,  # DAR0 from NONMEM
               A8=0, A9=0, A10=0, A11=0, A12=0,
               A13=150, A14=0, A15=0, A16=0)
    
    ode_func <- function(time, state, parms) {
      tgi_model(time, state, parms, dosing_times, input$dose_amt)
    }
    
    events <- list(
      data = data.frame(
        var = "A1",
        time = dosing_times,
        value = input$dose_amt * 1e6 / input$mw,
        method = "add"
      )
    )
    
    out <- as.data.frame(ode(
      y = state,
      times = times,
      func = ode_func,
      parms = parameters,
      events = events
    ))
    colnames(out)[2:(length(state)+1)] <- names(state)
    
    
    output$plotW <- renderPlot({
      plot(out$time, out$W, type = "l", col = "firebrick", lwd = 2,
           ylab = "Tumor Volume (mm³)", xlab = "Time")
    })
    output$plotADC <- renderPlot({
      plot(out$time, out$ADC, type = "l", col = "steelblue", lwd = 2,
           ylab = "ADC (ng/mL)", xlab = "Time")
    })
    output$plotPL <- renderPlot({
      plot(out$time, out$PLconc, type = "l", col = "forestgreen", lwd = 2,
           ylab = "PL (ng/mL)", xlab = "Time")
    })
    output$plotPLTotal <- renderPlot({
      plot(out$time, out$PL, type = "l", col = "purple", lwd = 2,
           ylab = "PL Total (nmol/kg)", xlab = "Time")
    })
    
    output$plotA8 <- renderPlot({
      plot(out$time, out$A8, type = "l", col = "darkgreen", lwd = 2,
           ylab = "A8 (tumor free ADC)", xlab = "Time")
    })
    output$plotA9 <- renderPlot({
      plot(out$time, out$A9, type = "l", col = "darkred", lwd = 2,
           ylab = "A9 (tumor bound ADC)", xlab = "Time")
    })
    output$plotA10 <- renderPlot({
      plot(out$time, out$A10, type = "l", col = "darkblue", lwd = 2,
           ylab = "A10 (PL tumor)", xlab = "Time")
    })
    output$plotA11 <- renderPlot({
      plot(out$time, out$A11, type = "l", col = "goldenrod", lwd = 2,
           ylab = "A11 (PL cell free)", xlab = "Time")
    })
    output$plotA12 <- renderPlot({
      plot(out$time, out$A12, type = "l", col = "brown", lwd = 2,
           ylab = "A12 (PL cell bound)", xlab = "Time")
    })
  
  output$plotAll <- renderPlot({
    par(mar=c(5,4,4,4))
    y1 <- out$W
    y1[y1 < 1e-4] <- 0
    y2 <- out$ADC
    y3 <- out$PLconc
    if (all(is.finite(y1)) && all(is.finite(y2)) && all(is.finite(y3))) {
      plot(out$time, y1, type = "l", col = "firebrick", lwd = 2,
           ylab = "Tumor Volume", xlab = "Time", ylim=c(0, max(y1)*1.2))
      par(new=TRUE)
      plot(out$time, y2, type = "l", col = "steelblue", lwd = 2, axes=FALSE,
           xlab="", ylab="", ylim=c(0, max(y2, y3)*1.2))
      lines(out$time, y3, col = "forestgreen", lwd = 2, lty=2)
      axis(side=4)
      mtext("PK Concentrations", side=4, line=3)
      legend("topright", legend=c("Tumor Volume", "ADC", "PL"),
             col=c("firebrick", "steelblue", "forestgreen"), lty=c(1,1,2), lwd=2)
    } else {
      plot.new()
      text(0.5, 0.5, "Invalid values detected (e.g., IC50 too small). Please adjust parameters.", cex = 1.2)
    }
  })
  })
}

shinyApp(ui = ui, server = server)
