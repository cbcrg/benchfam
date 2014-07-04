FROM fedora

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

RUN yum install -q -y which wget nano make gcc g++ gcc-gfortran expat-devel perl-CPAN perl-Net-SSLeay perl-IO-Socket-SSL openssl-devel unzip; \
  wget -q -O cpanm http://cpanmin.us; \
  chmod +x cpanm && mv cpanm bin/; \
  cpanm -q -n Env Net::SSLeay XML::Simple SOAP::Lite

RUN yum install -q -y python-setuptools
  
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
    rm satesrc-v2.2.7-2013Feb15.tar.gz; \
    ln -s /opt/satesrc-v2.2.7-2013Feb15 /opt/sate; 


RUN wget -q http://www.tcoffee.org/Packages/Beta/Version_10.00.04ad7ba/linux/T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64.tar.gz; \
  tar xf T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64.tar.gz; \
  mv T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64 /opt/tcoffee; \
  rm -rf T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64.tar.gz
  
ENV SATE_HOME /opt/sate/  
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/blast/bin:/opt/sate/sate-core/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/
ENV TEMP /tmp
ENV DIR_4_TCOFFEE /opt/tcoffee
ENV EMAIL_4_TCOFFEE tcoffee.msa@gmail.com
ENV CACHE_4_TCOFFEE /tmp/cache/
ENV LOCKDIR_4_TCOFFEE /tmp/lck/
ENV TMP_4_TCOFFEE /tmp/tmp/


