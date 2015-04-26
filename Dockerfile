FROM pditommaso/dkrbase:1.2

MAINTAINER Paolo Di Tommaso <paolo.ditommaso@gmail.com>

#
# Add Python setup tools 
#
RUN apt-get update --fix-missing && apt-get install -y python-setuptools
 
#
# BLAST
# 
RUN wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz -O- \
  | tar xz -C /opt && \
  ln -s /opt/ncbi-blast-2.2.29+/ /opt/blast; 

#
# MAFFT 
#  
RUN wget -q http://mafft.cbrc.jp/alignment/software/mafft-7.130-with-extensions-src.tgz -O- \
  | tar xz && \
  cd mafft-7.130-with-extensions/core && \
  make && make install && \
  cd $HOME && \
  rm -rf mafft-7.130-with-extensions 

#
# SATe
#
RUN easy_install -U dendropy && \
    cd /opt && \
    wget -q http://phylo.bio.ku.edu/software/sate/downloads2/src/satesrc-v2.2.7-2013Feb15.tar.gz -O- \
    | tar xz && \
    cd satesrc-v2.2.7-2013Feb15/sate-core/ && \
    python setup.py develop && \
    cd /opt && \
    ln -s /opt/satesrc-v2.2.7-2013Feb15 /opt/sate

#
# CLUSTALO
#
RUN wget -q http://www.clustal.org/omega/clustalo-1.2.0-Ubuntu-x86_64 && \
    chmod +x clustalo-1.2.0-Ubuntu-x86_64 && \
    mv clustalo-1.2.0-Ubuntu-x86_64 /usr/local/bin/clustalo

#
# T-COffee
#
RUN wget -q http://www.tcoffee.org/Packages/Beta/Version_10.00.04ad7ba/linux/T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64.tar.gz -O- \
  | tar xz && \
  mv T-COFFEE_installer_Version_10.00.04ad7ba_linux_x64 /opt/tcoffee 
  
#
# Perl modules required by T-Coffee
# 
RUN cpanm -q -n Env Net::SSLeay XML::Simple SOAP::Lite && \
  rm -rf /root/.cpanm/work/  
  
#
# SETUP the environment 
#
ENV SATE_HOME /opt/sate/  
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/blast/bin:/opt/sate/sate-core/bin:/opt/tcoffee/bin:/opt/tcoffee/plugins/linux/
ENV TEMP /tmp
ENV DIR_4_TCOFFEE /opt/tcoffee
ENV EMAIL_4_TCOFFEE tcoffee.msa@gmail.com
ENV CACHE_4_TCOFFEE /tmp/cache/
ENV LOCKDIR_4_TCOFFEE /tmp/lck/
ENV TMP_4_TCOFFEE /tmp/tmp/


