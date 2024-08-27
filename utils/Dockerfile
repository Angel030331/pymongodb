FROM mongo:latest

RUN useradd --create-home --shell /bin/bash user
WORKDIR /home/user
COPY requirements.txt ./

RUN apt-get update && apt-get install -y python3 python3-pip
RUN apt-get install -y libgtk-3-dev libwebkit2gtk-4.0-dev pkg-config

# ENV PKG_CONFIG_PATH=/usr/lib/x86_64-linux-gnu/pkgconfig

RUN pip install --no-cache-dir -r requirements.txt

USER user
COPY . /home/user/code

CMD ["tail", "-f", "/dev/null"]