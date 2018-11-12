envel.norm <- function(modelo=fit.model,iden=0,nome=seq(along = model.matrix(modelo)[,1]),sim=100,conf=.90,res=T,quad=T) {
  
  #
  # Descrição e detalhes:
  # A saída será o gráfico de probabilidade normal com envelopes simulados para um ajuste da distribuição normal.
  #
  # A opção res=F faz o gráfico de probabilidade meio-normal com envelopes simulados utilizando a distância de Cook,
  # possibilitando a detecção de pontos simultaneamente aberrantes e/ou influentes.
  #
  # Atenção: a função não funcionará corretamente se o ajuste possuir offsets! Neste caso é preciso adaptá-la como foi
  # feito na função envel.pois
  #
  # Os dados devem estar disponíveis pelo comando attach( ).
  #
  # Argumentos obrigatórios:
  # modelo: deve-se informar o objeto onde está o ajuste do modelo normal linear, caso não seja informado, a função
  # 	  procurará o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o número de observações que irá querer destacar. O padrão é não destacar ninguém (iden=0).
  #	Qualquer valor que não seja um inteiro positivo (por ex., negativo ou decimal) fará com que a função pergunte
  #	o número de pontos após a execução;
  # nome: esse argumento só é utilizado caso seja destacado algum ponto no gráfico. Caso não seja informado nada, os pontos
  #	identificados serão os números da ordem em que estão no banco de dados (os índices). Caso se queira, pode-se
  #	informar um vetor de nomes ou de identificações alternativas. Obrigatoriamente esse vetor deve ter o mesmo
  #	comprimento do banco de dados;
  # sim: número de simulações para gerar a banda de confiança. Atkinson sugere um mínimo de 20 simulações.
  #      O padrão é de 100;
  # conf: nível de confiança do envelope. O padrão é de 90%;
  # res: permite-se a escolha se o gráfico será feito com os resíduos (res=T, True, padrão) ou com a distância de Cook
  #      (res=F, False);
  # quad: o padrão (quad=T, True) faz um gráfico quadrado, enquanto quad=F (False) faz um gráfico utilizando a área máxima
  #       disponível.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4ª ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
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
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
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
    #Segundo McCullagh e Nelder (1989, pág.407) e Paula (2003, pág.57) deve-se usar qnorm((n+1:n+.5)/(2*n+1.125))
    #Segundo Neter et alli (1996, pág.597) deve-se usar qnorm((n+1:n-.125)/(2*n+0.5))
    qq<-qnorm((n+1:n+.5)/(2*n+1.125))
    plot(qq,sort(di),xlab="Quantil Meio-Normal",ylab="Distância de Cook", ylim=range(di,e1,e2), pch=16)
    nome<-nome[order(di)]
    r<-sort(di)
  } else {
    qq<-qnorm((1:n-.375)/(n+.25))
    plot(qq,sort(tsi),xlab="Quantil da Normal Padrão",ylab="Resíduo Padronizado", ylim=range(tsi,e1,e2), pch=16,main="Gráfico de envelope simulado com 95% de confiança")
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
  # Descrição e detalhes:
  # A saída terá seis gráficos:
  # 1º) Influência na Locação. O gráfico feito é das distâncias de Cook contra os valores ajustados. Utilizou-se o critério
  #     de destacar observações maiores do que duas vezes a média de todas as distâncias obtidas;
  # 2º) Influência Locação/Escala. A medida C, que é um aperfeiçoamento do DFFIT e também é conhecida como distância de Cook
  #     modificada, foi utilizada para medir a influência das observações nos parâmetros de locação e escala. O critério
  #     foi o de destacar observações maiores do que duas vezes a média de todas as distâncias obtidas;
  # 3º) Influência Local. A influência local consiste em procurar pontos que sob pequenas perturbações causam variações
  #     muito grandes nos resultados. O dmax é o autovetor que corresponde ao maior autovalor da matriz do processo de
  #     perturbações. Para maiores detalhes veja Paula (2003, págs.50-54 e 65-66). O critério foi o de destacar observações
  #     maiores do que duas vezes a média de todos os dmax's;
  # 4º) Pontos Alavanca. A matriz chapéu H=X%*%solve(t(X)%*%X)%*%t(X) é a matriz de projeção ortogonal de vetores no
  #     subespaço gerado pelas colunas da matrix X. Os pontos remotos nesse subespaço costumam ser considerados alavanca
  #     (leverage), por exercer uma forte influência no seu valor ajustado. Ou seja, esses pontos tem um perfil diferente
  #     dos demais com relação às variáveis explicativas. Ao fazer predições para um determinado vetor x, pode-se também
  #     obter a medida h para esse valor, através de h=t(x)%*%solve(t(X)%*%X)%*%x. Caso esse valor seja grande com relação
  #     aos pontos utilizados na estimação do modelo, isso é um indício em que a combinação de valores de x é uma
  #     extrapolação, mesmo que o valor separadamente de cada variável esteja dentro dos limites em que o modelo abrange;
  # 5º) Pontos Aberrantes. Um ponto é aberrante (discrepante, outlier) se o seu valor estiver mal ajustado pelo modelo.
  #     Como os resíduos utilizados tem uma distribuição t-Student de n-p-1 graus de liberdade, em que n é o número de
  #     observações e p o número de parâmetros, então adicionamos linhas tracejadas nos quantis 2.5% e 97.5%. Com isso
  #     esperamos que cerca de 5% dos pontos possam estar um pouco fora desses limites. Esse gráfico serve como indicação
  #     para detectar valores aberrantes marginalmente. Se o objetivo for detectar valores conjuntamente aberrantes
  #     deve-se construir o gráfico de envelopes ou utilizar critérios de comparações múltiplas, como o de Bonferroni
  #     que consistiria em utilizar os quantis 2.5%/n e 1-2.5%/n, uma vez que estamos fazendo n comparações. Para ter
  #     uma idéia, pode-se obter esses valores através dos comandos
  #	qt(.025/sum(summary(modelo)$df[1:2]),summary(modelo)$df[2]-1);
  #	qt(1-.025/sum(summary(modelo)$df[1:2]),summary(modelo)$df[2]-1);
  # 6º) Função de Variância. McCullagh e Nelder (1989, pág.400) sugere o gráfico dos resíduos absolutos contra os valores
  #     ajustados para checar se a função de variância adotada é adequada. O padrão esperado é de não encontrarmos nenhuma
  #     tendência. Funções de variância erradas irão resultar em tendências dos resíduos com a média. Tendências positivas
  #     indicam que a função de variância está crescendo muito devagar com a média, então deve-se aumentar a potência (no
  #     caso de uma função de variância da família potência). Uma linha suavizada pelo método lowess robusto é adicionada
  #     para ajudar na procura de tendências.
  #
  # Os dados devem estar disponíveis pelo comando attach( ).
  #
  # Argumentos obrigatórios:
  # modelo: deve-se informar o objeto onde está o ajuste do modelo normal linear, caso não seja informado, a
  # 	  função procurará o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o número de observações que irá querer destacar em cada gráfico. O vetor deve conter 6
  #	posições de números inteiros. A ordem que deve ser informada é a mesma em que os gráficos são feitos. Os
  #	componentes do vetor iguais a 0 indicam que não se quer que identifique pontos, se for um inteiro positivo irá
  #	automaticamente nos gráficos respectivos permitir que identifiquemos o número de pontos solicitados e qualquer
  #	outro valor (negativo ou decimal) parar nos gráficos e solicitar que especifiquemos o número de pontos a ser
  #	destacado. O padrão é c(0,0,0,0,0,0) caso não se entre com nada e c(-1,-1,-1,-1,-1,-1) caso se entre
  #	com qualquer coisa que não seja um vetor de 6 posições, como por ex.-1;
  # nome: esse argumento só é utilizado caso algum dos componentes do vetor da opção iden não seja 0. Caso não seja
  #	informado nada, os pontos identificados serão os números da ordem em que estão no banco de dados (índices).
  #	Caso se queira, pode-se informar um vetor de nomes ou de identificações alternativas. Obrigatoriamente
  #	esse vetor deve ter o mesmo comprimento do banco de dados.
  #
  # A função retorna os seguintes valores: ResPearsonStd, Di, Ci, Dmax e h.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referências:
  # MCCULLAGH, P. e NELDER, J. A. (1989). Generalized Linear Models. 2ª ed. Chapman and Hall, London.
  # NETER, J., KUTNER, M. H., NACHTSHEIM, C. J. and WASSERMAN, W. (1996). Applied Linear Statistical Models. 4ª ed.
  #    Mc Graw Hill, Boston.
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
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
  
  #para evitar divisão por 0 ao studentizar os residuos, mas tentando manter o valor exagerado da alavanca
  h[round(h,15)==1]<-0.999999999999999
  
  lms <- summary(modelo)
  s <- lms$sigma
  r <- resid(modelo)
  ts <- r/(s*sqrt(1-h))
  di <- (1/p)*(h/(1-h))*(ts^2)
  si <- lm.influence(modelo)$sigma
  tsi <- r/(si*sqrt(1-h))
  #dff <- sqrt(h/(1-h))*abs(tsi) #DFFIT
  ci <- sqrt( ((n-p)*h) / (p*(1-h)) )*abs(tsi) #aperfeiçoamento do DFFIT
  A <- diag(r)%*%H%*%diag(r)
  dmax <- abs(eigen(A)$vec[,1]/sqrt(eigen(A)$val[1]))
  m <- fitted(modelo)
  
  par(mfrow=c(2,3))
  
  plot(m,di,xlab="Valor Ajustado", ylab="Distância de Cook",main="Influência na Posição", ylim=c(0,max(di,2*mean(di))), pch=16)
  abline(2*mean(di),0,lty=2)
  while ( (!is.numeric(iden[1])) || (round(iden[1],0) != iden[1]) || (iden[1] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[1]<-as.numeric(out)
  }
  if(iden[1]>0) {identify(m,di,n=iden[1],labels=nome)}
  
  plot(m,ci,xlab="Valor Ajustado", ylab="Distância de Cook Modificada",main="Influência Posição/Escala", ylim=c(0,max(ci,2*mean(ci))), pch=16)
  abline(2*mean(ci),0,lty=2)
  while ( (!is.numeric(iden[2])) || (round(iden[2],0) != iden[2]) || (iden[2] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[2]<-as.numeric(out)
  }
  if(iden[2]>0) {identify(m,ci,n=iden[2],labels=nome)}
  
  plot(m,dmax,xlab="Valor Ajustado", ylab="dmax",main="Influência Local", ylim=c(0,max(dmax,2*mean(dmax))), pch=16)
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
  
  plot(m,tsi,xlab="Valor Ajustado", ylab="Resíduo Padronizado",main="Pontos Aberrantes", ylim=c(min(tsi)-1,max(tsi)+1), pch=16)
  abline(qt(.025,n-p-1),0,lty=2)
  abline(qt(1-.025,n-p-1),0,lty=2)
  while ( (!is.numeric(iden[5])) || (round(iden[5],0) != iden[5]) || (iden[5] < 0) ) {
    cat("Digite o num.de pontos a ser identificado (0=nenhum) e <enter> para continuar\n")
    out <- readline()
    iden[5]<-as.numeric(out)
  }
  if(iden[5]>0) {identify(m,tsi,n=iden[5],labels=nome)}
  
  plot(m,abs(tsi),xlab="Valor Ajustado", ylab="Resíduo Padronizado Absoluto",main="Função de Variância", pch=16)
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
  # Descrição e detalhes:
  # A saída terá gráficos do dmax para cada coeficiente de um modelo normal linear.
  #
  # A influência local consiste em procurar pontos que sob pequenas perturbações causam variações muito grandes nos
  # resultados. O dmax é o autovetor que corresponde ao maior autovalor da matriz do processo de perturbações. Para
  # maiores detalhes veja Paula (2003, págs.50-54). O critério foi o de destacar observações maiores do que duas
  # vezes a média de todos os dmax's absolutos.
  #
  # Os dados devem estar disponíveis pelo comando attach( ).
  #
  # Argumentos obrigatórios:
  # modelo: deve-se informar o objeto onde está o ajuste do modelo normal linear, caso não seja informado, a função
  # 	  procurará o ajuste no objeto fit.model;
  # 
  # Argumentos opcionais:
  # iden: caso deseje, informe o número de observações que irá querer destacar em cada gráfico. O vetor deve
  # 	conter números inteiros. A ordem que deve ser informada é a mesma das variáveis da opção var, caso seja
  #	utilizada, ou deve ter a mesma ordem da matriz modelo. Os componentes do vetor iguais a 0 indicam que não se
  #	quer que identifique pontos, se for um inteiro positivo irá automaticamente nos gráficos pertinentes permitir
  #	que identifiquemos o número de pontos solicitados e qualquer outro valor (negativo ou decimal) parará nos
  #	gráficos e irá solicitar que especifiquemos o número de pontos a ser destacado. O padrão é c(0,...,0) caso não
  #	se entre com nada e c(-1,...,-1) caso se entre com qualquer coisa que não satisfaça os requisitos citados
  #	de ser número inteiro, não negativo e de ter o mesmo comprimento da opção var ou da matriz modelo;
  # nome: esse argumento só é utilizado caso algum dos componentes do vetor da opção iden não seja 0. Caso não
  #	seja informado nada, os pontos identificados serão os números da ordem em que estão no banco de dados.
  #	Caso se queira, pode-se informar um vetor de nomes ou de identificações alternativas. Obrigatoriamente
  #	esse vetor deve ter o mesmo comprimento do banco de dados.
  #
  # Autor: Frederico Zanqueta Poleto <fred@poleto.com>, arquivo disponível em http://www.poleto.com
  #
  # Referência:
  # PAULA, G. A. (2003). Modelos de Regressão com apoio computacional. IME-USP, São Paulo. [Não publicado,
  #    disponível em http://www.ime.usp.br/~giapaula/Book.pdf]
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
    plot(dmax,xlab="Índice", ylab="dmax",main=dimnames(X)[[2]][i],ylim=c(min(-cut,dmax),max(cut,dmax)), pch=16)
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