default:
	@echo build3
	@echo run3
	@echo bash3
	@echo build4
	@echo run4
	@echo bash4
	@echo build5
	@echo run5
	@echo bash5
	@echo compose
	@echo rmd-look
	@echo purl

build3:
	docker build -q -t pshannon/tmsft3 --file Dockerfile.3 .

build4:
	docker build -q -t pshannon/tmsft4 --file Dockerfile.4 .

build5:
	docker build -q -t pshannon/tmsft5 --file Dockerfile.5 .



run3:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           -v /Users/paul/github/tms-brandLab/klf1-study/data:/home/rstudio/data \
           -v /Users/paul/github/tms-brandLab/klf1-study/code:/home/rstudio/code \
           -v /Users/paul/github/tms-brandLab/klf1-study/components:/home/rstudio/components \
           --name tmsft3 pshannon/tmsft3

run4:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           -e DISABLE_AUTH=true \
           -v /Users/paul/github/tms-brandLab/klf1-study/data:/home/rstudio/data \
           -v /Users/paul/github/tms-brandLab/klf1-study/code:/home/rstudio/code \
           -v /Users/paul/github/tms-brandLab/klf1-study/components:/home/rstudio/components \
           --name tmsft4 pshannon/tmsft4

run5:
	docker run --rm \
           -p 8888:8787 \
           -e PASSWORD=password \
           -e DISABLE_AUTH=true \
           -v /Users/paul/github/tms-brandLab/klf1-study/data:/home/rstudio/data \
           -v /Users/paul/github/tms-brandLab/klf1-study/code:/home/rstudio/code \
           -v /Users/paul/github/tms-brandLab/klf1-study/components:/home/rstudio/components \
           --name tmsft5 pshannon/tmsft5


bash3:
	docker exec -it tmsft3 bash

bash4:
	docker exec -it tmsft4 bash

bash5:
	docker exec -it tmsft5 bash

compose:
	docker-compose up

run-rstudio:
	open http://localhost:8889/

rmd-look:
	R -e "rmarkdown::render('klf1.Rmd')"
	open klf1.html

purl:
	R -e 'knitr::purl(input = "klf1.Rmd", output = "klf1.R",documentation = 0)'

