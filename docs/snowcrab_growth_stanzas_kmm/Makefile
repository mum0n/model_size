FILENAME=snowcrab_growth_stanzas_kmm
TEMPLATE=my_template.tex
PDFENGINE=lualatex
TMPDIR := $(shell mktemp -d)
MAKE_TMPDIR=TMDIR 

# export TEXINPUTS=.:media//:
# export BIBINPUTS=.:media//:
# export BSTINPUTS=.:media//:

all: latex pdf 

latex: 
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FILENAME).tex \
	--pdf-engine=$(PDFENGINE)

pdf: latex  
	pdflatex $(FILENAME).tex
	biber $(FILENAME) 
	pdflatex $(FILENAME).tex
	pdflatex $(FILENAME).tex

view: 
	zathura $(FILENAME).pdf 

test: 
	pandoc -s $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars+raw_tex \
	--to=latex --template=$(TEMPLATE) \
	--output=$(FILENAME).pdf \
	--pdf-engine=$(PDFENGINE)
	pdflatex $(FILENAME).tex

html:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=html5 \ 
	--output=$(FILENAME).html \
	--mathjax \
	--self-contained

epub:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=epub \
	--output=$(FILENAME).epub \
	--epub-cover-image=<cover-image> \
	--toc

docx:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=docx \
	--output=$(FILENAME).docx \


odt:
	pandoc $(FILENAME).md \
	--from=markdown+tex_math_single_backslash+tex_math_dollars \
	--to=odt \
	--output=$(FILENAME).odt \

  
git:
	git commit -m"update"
	git checkout master
	git merge develop
	git push
	git pull
	git checkout develop
	git merge master
	git status
	

 
watch: $(FILENAME).md 
	fswatch -o $^ | xargs -n1 -I{} make

.PHONY: clean all

clean:
	rm -rf *.aux *.bbl *.bcf *.blg *.log *.out  *.run.xml *.spl  *.docx *.odt *.epub *.html
	
