FROM python:3.8.12-bullseye
# RUN apt-get update
RUN useradd --create-home --shell /bin/bash app_user
WORKDIR /home/app_user
COPY requirements.txt /home/app_user/requirements.txt
RUN pip install --upgrade pip
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
USER app_user
COPY --chown=app_user:app_user . .
# COPY --chown=app_user:app_user ./mirnome-localizations-files ./for_localization_ref_files
CMD ["bash"]
