envel.norm <- function(modelo=fit.model,iden=0,nome=seq(along = model.matrix(modelo)[,1]),sim=100,conf=.90,res=T,quad=T) {
  
  #
  # Descri��o e detalhes:
  # A sa�da ser� o gr�fico de probabilidade normal com envelopes simulados para um ajuste da distribui��o normal.
  #
  # A op��o res=F faz o gr�fico de probabilidade meio-normal com envelopes simulados utilizando a dist�ncia de Cook,
  # possibilitando a detec��o de pontos simultaneamente aberrantes e/ou influentes.
  #
  # Aten��o: a fun��o n�o funcionar� corretamente se o ajuste possuir offsets! Neste caso � preciso adapt�-la como foi
  # feito na fun��o envel.pois
  #
  # Os dados devem estar dispon�veis pelo comando attach( ).
  #
  # Argumentos obrigat�rios:
  # modelo: deve-se informar o objeto onde est� o ajuste do modelo normal linear, caso n�o seja informado, a fun��o
  # 	  procurar� o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o n�mero de observa��es que ir� querer destacar. O padr�o � n�o destacar ningu�m (iden=0).
  #	Qualquer valor que n�o seja um inteiro positivo (por ex., negativo ou decimal) far� com que a fun��o pergunte
  #	o n�mero de pontos ap�s a execu��o;
  # nome: esse argumento s� � utilizado caso seja destacado algum ponto no gr�fico. Caso n�o seja informado nada, os pontos
  #	identificados ser�o os n�meros da ordem em que est�o no banco de dados (os �ndices). Caso se queira, pode-se
  #	informar um vetor de nomes ou de identifica��es alternativas. Obrigatoriamente esse vetor deve ter o mesmo
  #	comprimento do banco de dados;
  # sim: n�mero de simula��es para gerar a banda de confian�a. Atkinson sugere um m�nimo de 20 simula��es.
  #      O padr�o � de 100;
  # conf: n�vel de confian�a do envelope. O padr�o � de 90%;
  # res: permite-se a escolha se o gr�fico ser� feito com os res�duos (res=T, True, padr�o) ou com a dist�ncia de Cook
  #      (res=F, False);
  # quad: o padr�o (quad=T, True) faz um gr�fico quadrado, enquanto quad=F (False) faz um gr�fico utilizando a �rea m�xima
  #       dispon�vel.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo dispon�vel em http://www.poleto.com
  #
  # Refer�ncias:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2� ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4� ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regress�o com apoio computacional. IME-USP, S�o Paulo. [N�o publicado,
  #    dispon�vel em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  # Exemplos:
  # envel.norm(ajuste,sim=10000,conf=.95)
  # envel.norm(ajuste,res=F)
  #
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  alfa<-(1-conf)/2
  X <- model.matrix(modelo)
  y<-predict(modelo)+resid(modelo)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  m <- fitted(modelo)
  
  #para evitar divis�o por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  si <- lm.influence(modelo)$sigma
  r <- resid(modelo)
  tsi <- r/(si*sqrt(1-h))
  sigma<-summary(modelo)$sigma
  ti <- r/(sigma*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ti^2)
  
  e <- matrix(0,n,sim)
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:sim) {
    resp <- rnorm(n,m,sigma)
    fit <- lm(resp~X-1)
    ti<-resid(fit)/(summary(fit)$sigma*sqrt(1-h))
    if(res==F) {
      e[,i] <- (1/p)*(h/(1-h))*(ti^2)
    } else {
      e[,i] <- ti*sqrt( (n-p-1)/(n-p-(ti^2)) )
    }	
    e[,i] <- sort(e[,i])
  }
  
  for(i in 1:n) {
    eo <- sort(e[i,])
    e1[i] <- quantile(eo,alfa)
    e2[i] <- quantile(eo,1-alfa)
  }
  
  med <- apply(e,1,median)
  
  if(quad==T) {
    par(pty="s")
  }
  if(res==F) {
    #Segundo McCullagh e Nelder (1989, p�g.407) e Paula (2003, p�g.57) deve-se usar qnorm((n+1:n+.5)/(2*n+1.125))
    #Segundo Neter et alli (1996, p�g.597) deve-se usar qnorm((n+1:n-.125)/(2*n+0.5))
    qq<-qnorm((n+1:n+.5)/(2*n+1.125))
    plot(qq,sort(di),xlab="Quantil Meio-Normal",ylab="Dist�ncia de Cook", ylim=range(di,e1,e2), pch=16)
    nome<-nome[order(di)]
    r<-sort(di)
  } else {
    qq<-qnorm((1:n-.375)/(n+.25))
    plot(qq,sort(tsi),xlab="Quantil da Normal Padr�o",ylab="Res�duo Padronizado", ylim=range(tsi,e1,e2), pch=16,main="Gr�fico de envelope simulado com 95% de confian�a")
    nome<-nome[order(tsi)]
    r<-sort(tsi)
  }
  lines(qq,e1,lty=1)
  lines(qq,e2,lty=1)
  lines(qq,med,lty=2)
  while ( (!is.numeric(iden)) || (round(iden,0) != iden) || (iden < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden<-as.numeric(out)
  }
  if(iden>0) {identify(qq,r,n=iden,labels=nome)}
  if(quad==T) {
    par(pty="m")
  }
  cat("\nBanda de ",conf*100,"% de confianca, obtida por ",sim," simulacoes.\n")
}
diag.norm <- function(modelo=fit.model,iden=c(0,0,0,0,0,0),nome=seq(along = model.matrix(modelo)[,1])) {
  
  #
  # Descri��o e detalhes:
  # A sa�da ter� seis gr�ficos:
  # 1�) Influ�ncia na Loca��o. O gr�fico feito � das dist�ncias de Cook contra os valores ajustados. Utilizou-se o crit�rio
  #     de destacar observa��es maiores do que duas vezes a m�dia de todas as dist�ncias obtidas;
  # 2�) Influ�ncia Loca��o/Escala. A medida C, que � um aperfei�oamento do DFFIT e tamb�m � conhecida como dist�ncia de Cook
  #     modificada, foi utilizada para medir a influ�ncia das observa��es nos par�metros de loca��o e escala. O crit�rio
  #     foi o de destacar observa��es maiores do que duas vezes a m�dia de todas as dist�ncias obtidas;
  # 3�) Influ�ncia Local. A influ�ncia local consiste em procurar pontos que sob pequenas perturba��es causam varia��es
  #     muito grandes nos resultados. O dmax � o autovetor que corresponde ao maior autovalor da matriz do processo de
  #     perturba��es. Para maiores detalhes veja Paula (2003, p�gs.50-54 e 65-66). O crit�rio foi o de destacar observa��es
  #     maiores do que duas vezes a m�dia de todos os dmax's;
  # 4�) Pontos Alavanca. A matriz chap�u H=X%*%solve(t(X)%*%X)%*%t(X) � a matriz de proje��o ortogonal de vetores no
  #     subespa�o gerado pelas colunas da matrix X. Os pontos remotos nesse subespa�o costumam ser considerados alavanca
  #     (leverage), por exercer uma forte influ�ncia no seu valor ajustado. Ou seja, esses pontos tem um perfil diferente
  #     dos demais com rela��o �s vari�veis explicativas. Ao fazer predi��es para um determinado vetor x, pode-se tamb�m
  #     obter a medida h para esse valor, atrav�s de h=t(x)%*%solve(t(X)%*%X)%*%x. Caso esse valor seja grande com rela��o
  #     aos pontos utilizados na estima��o do modelo, isso � um ind�cio em que a combina��o de valores de x � uma
  #     extrapola��o, mesmo que o valor separadamente de cada vari�vel esteja dentro dos limites em que o modelo abrange;
  # 5�) Pontos Aberrantes. Um ponto � aberrante (discrepante, outlier) se o seu valor estiver mal ajustado pelo modelo.
  #     Como os res�duos utilizados tem uma distribui��o t-Student de n-p-1 graus de liberdade, em que n � o n�mero de
  #     observa��es e p o n�mero de par�metros, ent�o adicionamos linhas tracejadas nos quantis 2.5% e 97.5%. Com isso
  #     esperamos que cerca de 5% dos pontos possam estar um pouco fora desses limites. Esse gr�fico serve como indica��o
  #     para detectar valores aberrantes marginalmente. Se o objetivo for detectar valores conjuntamente aberrantes
  #     deve-se construir o gr�fico de envelopes ou utilizar crit�rios de compara��es m�ltiplas, como o de Bonferroni
  #     que consistiria em utilizar os quantis 2.5%/n e 1-2.5%/n, uma vez que estamos fazendo n compara��es. Para ter
  #     uma id�ia, pode-se obter esses valores atrav�s dos comandos
  #	qt(.025/sum(summary(modelo)$df[1:2]),summary(modelo)$df[2]-1);
  #	qt(1-.025/sum(summary(modelo)$df[1:2]),summary(modelo)$df[2]-1);
  # 6�) Fun��o de Vari�ncia. McCullagh e Nelder (1989, p�g.400) sugere o gr�fico dos res�duos absolutos contra os valores
  #     ajustados para checar se a fun��o de vari�ncia adotada � adequada. O padr�o esperado � de n�o encontrarmos nenhuma
  #     tend�ncia. Fun��es de vari�ncia erradas ir�o resultar em tend�ncias dos res�duos com a m�dia. Tend�ncias positivas
  #     indicam que a fun��o de vari�ncia est� crescendo muito devagar com a m�dia, ent�o deve-se aumentar a pot�ncia (no
  #     caso de uma fun��o de vari�ncia da fam�lia pot�ncia). Uma linha suavizada pelo m�todo lowess robusto � adicionada
  #     para ajudar na procura de tend�ncias.
  #
  # Os dados devem estar dispon�veis pelo comando attach( ).
  #
  # Argumentos obrigat�rios:
  # modelo: deve-se informar o objeto onde est� o ajuste do modelo normal linear, caso n�o seja informado, a
  # 	  fun��o procurar� o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o n�mero de observa��es que ir� querer destacar em cada gr�fico. O vetor deve conter 6
  #	posi��es de n�meros inteiros. A ordem que deve ser informada � a mesma em que os gr�ficos s�o feitos. Os
  #	componentes do vetor iguais a 0 indicam que n�o se quer que identifique pontos, se for um inteiro positivo ir�
  #	automaticamente nos gr�ficos respectivos permitir que identifiquemos o n�mero de pontos solicitados e qualquer
  #	outro valor (negativo ou decimal) parar nos gr�ficos e solicitar que especifiquemos o n�mero de pontos a ser
  #	destacado. O padr�o � c(0,0,0,0,0,0) caso n�o se entre com nada e c(-1,-1,-1,-1,-1,-1) caso se entre
  #	com qualquer coisa que n�o seja um vetor de 6 posi��es, como por ex.-1;
  # nome: esse argumento s� � utilizado caso algum dos componentes do vetor da op��o iden n�o seja 0. Caso n�o seja
  #	informado nada, os pontos identificados ser�o os n�meros da ordem em que est�o no banco de dados (�ndices).
  #	Caso se queira, pode-se informar um vetor de nomes ou de identifica��es alternativas. Obrigatoriamente
  #	esse vetor deve ter o mesmo comprimento do banco de dados.
  #
  # A fun��o retorna os seguintes valores: ResPearsonStd, Di, Ci, Dmax e h.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo dispon�vel em http://www.poleto.com
  #
  # Refer�ncias:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2� ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4� ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regress�o com apoio computacional. IME-USP, S�o Paulo. [N�o publicado,
  #    dispon�vel em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  # Exemplos:
  # diag.norm(ajuste,iden=c(1,5,2,4,3,0),nome=estados)
  # diag.norm(ajuste,iden=-1)
  #
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  if(length(iden)<6) {
    iden<-c(-1,-1,-1,-1,-1,-1)
  }
  
  X <- model.matrix(modelo)
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%solve(t(X)%*%X)%*%t(X)
  h <- diag(H)
  
  #para evitar divis�o por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  lms <- summary(modelo)
  s <- lms$sigma
  r <- resid(modelo)
  ts <- r/(s*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ts^2)
  si <- lm.influence(modelo)$sigma
  tsi <- r/(si*sqrt(1-h))
  #dff <- sqrt(h/(1-h))*abs(tsi) #DFFIT
  ci <- sqrt( ((n-p)*h) / (p*(1-h)) )*abs(tsi) #aperfei�oamento do DFFIT
  A <- diag(r)%*%H%*%diag(r)
  dmax <- abs(eigen(A)$vec[,1]/sqrt(eigen(A)$val[1]))
  m <- fitted(modelo)
  
  par(mfrow=c(2,3))
  
  plot(m,di,xlab="Valor Ajustado", ylab="Dist�ncia de Cook",main="Influ�ncia na Posi��o", ylim=c(0,max(di,2*mean(di))), pch=16)
  abline(2*mean(di),0,lty=2)
  while ( (!is.numeric(iden[1])) || (round(iden[1],0) != iden[1]) || (iden[1] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[1]<-as.numeric(out)
  }
  if(iden[1]>0) {identify(m,di,n=iden[1],labels=nome)}
  
  plot(m,ci,xlab="Valor Ajustado", ylab="Dist�ncia de Cook Modificada",main="Influ�ncia Posi��o/Escala", ylim=c(0,max(ci,2*mean(ci))), pch=16)
  abline(2*mean(ci),0,lty=2)
  while ( (!is.numeric(iden[2])) || (round(iden[2],0) != iden[2]) || (iden[2] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[2]<-as.numeric(out)
  }
  if(iden[2]>0) {identify(m,ci,n=iden[2],labels=nome)}
  
  plot(m,dmax,xlab="Valor Ajustado", ylab="dmax",main="Influ�ncia Local", ylim=c(0,max(dmax,2*mean(dmax))), pch=16)
  abline(2*mean(dmax),0,lty=2)
  while ( (!is.numeric(iden[3])) || (round(iden[3],0) != iden[3]) || (iden[3] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[3]<-as.numeric(out)
  }
  if(iden[3]>0) {identify(m,dmax,n=iden[3],labels=nome)}
  
  plot(m,h,xlab="Valor Ajustado", ylab="Medida h",main="Pontos Alavanca", ylim=c(0,max(h,2*p/n)), pch=16)
  abline(2*p/n,0,lty=2)
  while ( (!is.numeric(iden[4])) || (round(iden[4],0) != iden[4]) || (iden[4] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[4]<-as.numeric(out)
  }
  if(iden[4]>0) {identify(m,h,n=iden[4],labels=nome)}
  
  plot(m,tsi,xlab="Valor Ajustado", ylab="Res�duo Padronizado",main="Pontos Aberrantes", ylim=c(min(tsi)-1,max(tsi)+1), pch=16)
  abline(qt(.025,n-p-1),0,lty=2)
  abline(qt(1-.025,n-p-1),0,lty=2)
  while ( (!is.numeric(iden[5])) || (round(iden[5],0) != iden[5]) || (iden[5] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[5]<-as.numeric(out)
  }
  if(iden[5]>0) {identify(m,tsi,n=iden[5],labels=nome)}
  
  plot(m,abs(tsi),xlab="Valor Ajustado", ylab="Res�duo Padronizado Absoluto",main="Fun��o de Vari�ncia", pch=16)
  lines(lowess(m,abs(tsi)))
  while ( (!is.numeric(iden[6])) || (round(iden[6],0) != iden[6]) || (iden[6] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[6]<-as.numeric(out)
  }
  if(iden[6]>0) {identify(m,abs(tsi),n=iden[6],labels=nome)}
  
  par(mfrow=c(1,1))
  list(ResPearsonStd=tsi,Di=di,Ci=ci,Dmax=dmax,h=h)
}
dmax.norm <- function(modelo=fit.model,iden=double(ncol(model.matrix(modelo))),nome=seq(along = model.matrix(modelo)[,1])) {
  
  #
  # Descri��o e detalhes:
  # A sa�da ter� gr�ficos do dmax para cada coeficiente de um modelo normal linear.
  #
  # A influ�ncia local consiste em procurar pontos que sob pequenas perturba��es causam varia��es muito grandes nos
  # resultados. O dmax � o autovetor que corresponde ao maior autovalor da matriz do processo de perturba��es. Para
  # maiores detalhes veja Paula (2003, p�gs.50-54). O crit�rio foi o de destacar observa��es maiores do que duas
  # vezes a m�dia de todos os dmax's absolutos.
  #
  # Os dados devem estar dispon�veis pelo comando attach( ).
  #
  # Argumentos obrigat�rios:
  # modelo: deve-se informar o objeto onde est� o ajuste do modelo normal linear, caso n�o seja informado, a fun��o
  # 	  procurar� o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o n�mero de observa��es que ir� querer destacar em cada gr�fico. O vetor deve
  # 	conter n�meros inteiros. A ordem que deve ser informada � a mesma das vari�veis da op��o var, caso seja
  #	utilizada, ou deve ter a mesma ordem da matriz modelo. Os componentes do vetor iguais a 0 indicam que n�o se
  #	quer que identifique pontos, se for um inteiro positivo ir� automaticamente nos gr�ficos pertinentes permitir
  #	que identifiquemos o n�mero de pontos solicitados e qualquer outro valor (negativo ou decimal) parar� nos
  #	gr�ficos e ir� solicitar que especifiquemos o n�mero de pontos a ser destacado. O padr�o � c(0,...,0) caso n�o
  #	se entre com nada e c(-1,...,-1) caso se entre com qualquer coisa que n�o satisfa�a os requisitos citados
  #	de ser n�mero inteiro, n�o negativo e de ter o mesmo comprimento da op��o var ou da matriz modelo;
  # nome: esse argumento s� � utilizado caso algum dos componentes do vetor da op��o iden n�o seja 0. Caso n�o
  #	seja informado nada, os pontos identificados ser�o os n�meros da ordem em que est�o no banco de dados.
  #	Caso se queira, pode-se informar um vetor de nomes ou de identifica��es alternativas. Obrigatoriamente
  #	esse vetor deve ter o mesmo comprimento do banco de dados.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo dispon�vel em http://www.poleto.com
  #
  # Refer�ncia:
  # PAULA, G. A. (2003). Modelos de Regress�o com apoio computacional. IME-USP, S�o Paulo. [N�o publicado,
  #    dispon�vel em http://www.ime.usp.br/~giapaula/Book.pdf]
  #
  # Exemplo:
  # dmax.norm(ajuste,iden=-1,nome=estados)
  #
  
  if( class(modelo)[1]=="lm" || (class(modelo)[1]=="glm" && (modelo$family[[1]]=="Gaussian" | modelo$family[[1]]=="gaussian")) ) {
    
  } else {
    stop(paste("\nA classe do objeto deveria ser lm ou glm (com distribuicao gaussian) !!!"))
  }
  
  X <- model.matrix(modelo)
  n <- nrow(X)
  p <- ncol(X)
  r<-resid(modelo)
  
  if(p>length(iden)) {
    iden<-rep(-1,p)
  }
  
  if (p>2) {
    if (p>8) {	
      par(mfrow=c(3,ceiling(p/3)))
    } else {
      par(mfrow=c(2,ceiling(p/2)))
    }
  } else {
    par(mfrow=c(1,ceiling(p)))
  }
  
  for(i in 1:p) {
    expl<-""
    expli<-0
    for(j in 1:p) {
      if(j!=i) {
        if(expli>0) {
          expli<-expli+1
          expl<-paste(expl,"+X[,",j,"]",sep="")
        } else {
          expli<-1
          expl<-paste("X[,",j,"]",sep="")
        }
      }
    }
    v<-resid(lm(as.formula(paste("X[,",i,"]~",expl,"-1"))))
    dmax<-v*r
    dmax<-dmax/sqrt( t(dmax)%*%dmax )
    #cut<-2/sqrt(n)
    cut<-2*mean(abs(dmax))
    plot(dmax,xlab="�ndice", ylab="dmax",main=dimnames(X)[[2]][i],ylim=c(min(-cut,dmax),max(cut,dmax)), pch=16)
    abline(-cut,0,lty=2)
    abline(cut,0,lty=2)
    while ( (!is.numeric(iden[i])) || (round(iden[i],0) != iden[i]) || (iden[i] < 0) ) {
      cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
      out <- readline()
      iden[i]<-as.numeric(out)
    }
    if(iden[i]>0) {identify(dmax,n=iden[i],labels=nome)}
  }
  
  par(mfrow=c(1,1))
  cat("\n")
}