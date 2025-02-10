default:
	@echo build0
	@echo build1
	@echo run
	@echo stop
	@echo bash
	@echo compose

build:
	docker build -t pshannon/tms --file Dockerfile .

build0:
	docker build -t tms0 --file Dockerfile0 .

build1:
	docker build -t tms1 --file Dockerfile1 .


run:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           --name tms tms1

rstudio:
	docker run --rm -ti -e PASSWORD=tms -e ROOT=true -p 8787:8787 tms1

compose:
	docker-compose up

bash:
	docker exec -it tms1 bash

stop:
	docker rm -f tms1
