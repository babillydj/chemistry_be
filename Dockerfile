# pull official base image
FROM python:slim

# set work directory
WORKDIR /src/app

# set environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# install built-in dependency
RUN apt-get update && apt-get -y dist-upgrade
RUN apt-get -y install build-essential libssl-dev libffi-dev libblas3 libc6 liblapack3 gcc python3-dev python3-pip cython3
RUN apt-get -y install python3-numpy python3-scipy 

# install miniconda for redkit
RUN apt install netcat-traditional
RUN apt-get install -y --no-install-recommends \
    git \
    wget \
    g++ \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN apt-get install netcat
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh \
    && echo "Running $(conda --version)"

# install dependencies
RUN pip install --upgrade pip
RUN pip install poetry
COPY pyproject.toml /src/app/
RUN poetry config virtualenvs.create false \
  && poetry install --no-interaction --no-ansi


# copy entrypoint.sh
COPY entrypoint.sh .
RUN sed -i 's/\r$//g' /src/app/entrypoint.sh
RUN chmod +x /src/app/entrypoint.sh

# copy project
COPY . .

# run entrypoint.sh
ENTRYPOINT ["/src/app/entrypoint.sh"]
