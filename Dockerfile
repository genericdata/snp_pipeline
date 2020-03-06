FROM centos:centos7.7.1908

RUN yum -y install \
	git \
	wget \
	java-1.8.0-openjdk \
	java-1.8.0-openjdk-devel

ENV APPS_ROOT /apps
RUN mkdir -p ${APPS_ROOT}

###############################################
#BWA = 'bwa/intel/0.7.17'

ENV BWA_VERSION 0.7.17

#ENV BWA_ROOT ${APPS_ROOT}/bwa/${BWA_VERSION}
ENV BWA_HOME ${APPS_ROOT}/bwa/${BWA_VERSION}
ENV PATH ${BWA_HOME}/bin:${PATH}

RUN git clone --branch v${BWA_VERSION} https://github.com/lh3/bwa.git ${BWA_HOME}
RUN cd ${BWA_HOME}; make

###############################################
#PICARD = 'picard/2.17.11'
# requires java 1.8
#RUN yum install -y \
 #      java-1.8.0-openjdk \
 #      java-1.8.0-openjdk-devel

ENV PICARD_VERSION 2.17.11

ENV JAVA_HOME /etc/alternatives/jre
#ENV PICARD_ROOT ${APPS_ROOT}/picard/${PICARD_VERSION}
ENV PICARD_HOME ${APPS_ROOT}/picard/${PICARD_VERSION}
ENV PICARD_JAR ${PICARD_HOME}/picard-${PICARD_VERSION}.jar

RUN wget https://github.com/broadinstitute/picard/releases/download/2.17.11/picard.jar -O ${PICARD_HOME}

#GATK = 'gatk/4.1.3.0'
###############################################
#R = 'r/intel/3.4.2'
###############################################
#SAMTOOLS = 'samtools/intel/1.9'
###############################################
#SNPEFF = 'snpeff/4.3i'


