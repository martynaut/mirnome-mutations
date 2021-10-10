FROM python:3.8.12-bullseye
# FROM conda/c3i-linux-64:latest
# RUN apt-get update
RUN useradd --create-home --shell /bin/bash app_user
WORKDIR /home/app_user
# RUN conda create --name mirmut_env  python=3.8
# SHELL ["conda", "run", "-n", "mirmut_env", "/bin/bash", "-c"]
COPY requirements.txt /home/app_user/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt
# WORKDIR /home
# RUN wget https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-2.4.18.tar.gz
# RUN tar -zxvf ViennaRNA-2.4.18.tar.gz
# WORKDIR /home/ViennaRNA-2.4.18
# RUN ./configure
# RUN make
# RUN make check
# RUN make install
RUN pip install ViennaRNA
WORKDIR /home/app_user
# RUN chown app_user:app_user /home/app_user/.bashrc
# SHELL ["/bin/bash", "-c"]
# RUN export PATH=$/home/ViennaRNA-2.4.18:${PATH}
# RUN source /home/app_user/.bashrc
USER app_user
COPY --chown=app_user:app_user . .
COPY --chown=app_user:app_user ./mirnome-localizations-files ./for_localization_ref_files
CMD ["bash"]
