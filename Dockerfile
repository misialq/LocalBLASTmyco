FROM python:3.6.5
LABEL maintainer "ziomal.mike@gmail.com"

RUN apt-get update && apt-get install -y wget

RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
RUN tar xzf ncbi-blast-2.9.0+-x64-linux.tar.gz

ENV PATH=".:/ncbi-blast-2.9.0+/bin:${PATH}"

COPY Pipfile /Pipfile
COPY Pipfile.lock /Pipfile.lock

RUN pip install --upgrade pip
RUN pip install pipenv
RUN pipenv install --system --deploy

COPY LocalBLAST.py /LocalBLAST.py
COPY LocalBLAST/ /LocalBLAST/
COPY blast/ /blast/

ENTRYPOINT ["python","LocalBLAST.py"]