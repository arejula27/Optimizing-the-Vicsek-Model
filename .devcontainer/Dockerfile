FROM ubuntu:latest

RUN apt update -y
RUN apt upgrade -y

RUN apt install git -y
RUN apt install g++ -y
RUN apt install make -y
RUN apt install zsh -y
RUN apt install sudo -y
RUN apt install clang-format -y
RUN apt install gnuplot -y


RUN echo '%sudo ALL=(ALL) NOPASSWD:ALL' >>/etc/sudoers
RUN useradd --create-home --password vscode --shell /bin/zsh vscode
COPY ./.zshrc /home/vscode


