FROM mongo:latest

RUN useradd --create-home --shell /bin/bash/user
WORKDIR /home/user
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

USER user
COPY . /home/user/code

CMD ["bash"]