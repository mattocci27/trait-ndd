FROM ubuntu:20.04

RUN  sed -i.bak -e "s%http://archive.ubuntu.com/%http://mirrors.tuna.tsinghua.edu.cn/%g" /etc/apt/sources.list \
  && sed -i.bak -e "s%http://security.ubuntu.com/%http://mirrors.tuna.tsinghua.edu.cn/%g" /etc/apt/sources.list

ENV TZ=Asia/Shanghai

RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


ENV R_VERSION=${R_VERSION:-3.6.3} \
    CRAN=${CRAN:-https//mirror.lzu.edu.cn/CRAN} \ 
    LC_ALL=en_US.UTF-8 \
    LANG=en_US.UTF-8 \
    TERM=xterm

RUN apt-get update && apt-get -y install lsb-release locales

RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
  && locale-gen en_US.utf8 \
  && /usr/sbin/update-locale LANG=en_US.UTF-8

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    bash-completion \
    ca-certificates \
    devscripts \
    file \
    fonts-texgyre \
    g++ \
    gfortran \
    gsfonts \
    libblas-dev \
    libbz2-* \
    libcurl4 \
    libicu* \
    libpcre2* \
    libjpeg-turbo* \
    libopenblas-openmp-dev \
    libpangocairo-* \
    libpng16* \
    libreadline8 \
    libtiff* \
    liblzma* \
    locales \
    make \
    unzip \
    zip \
    zlib1g


RUN apt-get update && apt-get install -y curl \
    default-jdk \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libicu-dev \
    libpcre2-dev \
    libpng-dev \
    libreadline-dev \
    libtiff5-dev \
    liblzma-dev \
    libx11-dev \
    libxt-dev \
    perl \
    rsync \
    subversion \
    tcl-dev \
    tk-dev \
    texinfo \
    texlive-extra-utils \
    texlive-fonts-recommended \
    texlive-fonts-extra \
    texlive-latex-recommended \
    texlive-latex-extra \
    x11proto-core-dev \
    xauth \
    xfonts-base \
    xvfb \
    wget \
    zlib1g-dev

RUN apt-get update && apt-get install -y xpra xserver-xephyr xinit xauth xclip x11-xserver-utils
