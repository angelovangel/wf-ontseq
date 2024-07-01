FROM aangeloo/nxf-tgs:latest

LABEL author="Angel Angelov <aangeloo@gmail.com>"
LABEL description="Docker image for the wf-ontseq script"

COPY ./ont-plasmid.sh .
RUN chmod +x ont-plasmid.sh

ENTRYPOINT [ "/ont-plasmid.sh" ]

# usa as:
# 
