default:
	@echo build3
	@echo run3
	@echo bash3
	@echo build

build3:
	docker build -q -t pshannon/tmsft3 --file Dockerfile.3 .

build:
	docker build -q -t pshannon/tmsft --file Dockerfile .



run3:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           -v /Users/paul/github/tms-brandLab/klf1-study/data:/home/rstudio/data \
           -v /Users/paul/github/tms-brandLab/klf1-study/code:/home/rstudio/code \
           -v /Users/paul/github/tms-brandLab/klf1-study/components:/home/rstudio/components \
           --name tmsft3 pshannon/tmsft3

run:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           -e DISABLE_AUTH=true \
           -v /Users/paul/github/tms-brandLab/klf1-study/data:/home/rstudio/data \
           -v /Users/paul/github/tms-brandLab/klf1-study/code:/home/rstudio/code \
           -v /Users/paul/github/tms-brandLab/klf1-study/components:/home/rstudio/components \
           --name tmsft pshannon/tmsft


bash3:
	docker exec -it tmsft3 bash

bash:
	docker exec -it tmsft bash

