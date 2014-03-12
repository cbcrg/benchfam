FROM fedora

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN yum install -q -y wget nano make gcc g++ gcc-gfortran expat-devel perl-CPAN perl-Net-SSLeay perl-IO-Socket-SSL openssl-devel unzip python-setuptools; \
  wget -q -O cpanm http://cpanmin.us; \
  chmod +x cpanm && mv cpanm bin/; \
  cpanm -q -n Net::SSLeay XML::Simple SOAP::Lite
  
RUN wget -q --no-cookies --no-check-certificate --header "Cookie: gpw_e24=http%3A%2F%2Fwww.oracle.com%2F" \
  http://download.oracle.com/otn-pub/java/jdk/7u51-b13/jdk-7u51-linux-x64.rpm -O jdk-7u51-linux-x64.rpm; \
  rpm -vhi jdk-7u51-linux-x64.rpm; \
  rm jdk-*
    
RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.29+-x64-linux.tar.gz; \
  tar xf ncbi-blast-2.2.29+-x64-linux.tar.gz; \
  mv ncbi-blast-2.2.29+ /opt/; \
  rm -rf ncbi-*; \
  ln -s /opt/ncbi-blast-2.2.29+/ /opt/blast; 
  
RUN wget -q http://mafft.cbrc.jp/alignment/software/mafft-7.130-with-extensions-src.tgz; \
  tar xf mafft-7.130-with-extensions-src.tgz; \
  cd mafft-7.130-with-extensions/core; \
  make;  \
  make install; \
  cd $HOME; \
  rm -rf mafft-7.130-with-extensions; \
  rm mafft-7.130-with-extensions-src.tgz; 

RUN easy_install -U dendropy; \
    cd /opt; \
    wget -q http://phylo.bio.ku.edu/software/sate/downloads2/src/satesrc-v2.2.7-2013Feb15.tar.gz; \
    tar xf satesrc-v2.2.7-2013Feb15.tar.gz; \
    cd satesrc-v2.2.7-2013Feb15/sate-core/; \
    python setup.py develop; \
    cd /opt; \
    ln -s /opt/satesrc-v2.2.7-2013Feb15/sate-core/ /opt/sate; 

RUN wget -q http://www.tcoffee.org/Packages/Stable/Latest/linux/T-COFFEE_installer_Version_10.00.r1613_linux_x64.bin; \
  chmod +x T-COFFEE_*; \
  ./T-COFFEE_installer_Version_10.00.r1613_linux_x64.bin --mode unattended --user_email tcoffee.msa@gmail.com --installdir /opt/tcoffee; \
  rm -rf T-COFFEE_*; \
  rm -rf .bash*; 
  
  
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/blast/bin:/opt/sate/bin:/opt/tcoffee/bin
ENV DIR_4_TCOFFEE /opt/tcoffee
ENV CACHE_4_TCOFFEE /.t_coffee/cache/
ENV MAFFT_BINARIES /opt/tcoffee/plugins/linux/
ENV EMAIL_4_TCOFFEE tcoffee.msa@gmail.com
ENV LOCKDIR_4_TCOFFEE /opt/tcoffee/lck/
ENV TMP_4_TCOFFEE /opt/tcoffee/tmp/