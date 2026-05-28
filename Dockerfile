# Pin digest so builds work offline when the base image is already pulled.
# If Hub auth fails, run: docker pull rocker/r-ver:4.4.1 && docker build --pull=false -t tf-heatmap-shiny .
FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libsodium-dev \
    libsass-dev \
    libuv1-dev \
    cmake \
    g++ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY install.R /app/install.R
RUN Rscript install.R

COPY app.R /app/app.R
COPY R/ /app/R/
COPY test/ /app/test/
COPY GSAD_DAPseq_TFs.R /app/GSAD_DAPseq_TFs.R

EXPOSE 3838

ENV SHINY_PORT=3838

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=as.integer(Sys.getenv('SHINY_PORT', '3838')))"]
