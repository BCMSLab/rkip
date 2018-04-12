R_OPT=--vanilla

SCRIPTS=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

DATA=data
LOG=log

all: $(DATA)/*.tsv $(DATA)/rkip_wgcna.rda $(FIG_DIR)/* $(TAB_DIR)/* clean

# install/update the functions
install:
	Rscript -e "devtools::install('..')"

# make directory tree
dir_tree:
	test ! -d $(LOG) && mkdir $(LOG) || echo '$(LOG)/ is already there.'
	test ! -d $(DATA) && mkdir $(DATA) || echo '$(DATA)/ is already there.'
	test ! -d $(MANUSCRIPT) && mkdir $(MANUSCRIPT) || echo '$(MANUSCRIPT)/ is already there.'
	test ! -d $(TAB_DIR) && mkdir $(TAB_DIR) || echo '$(TAB_DIR)/ is already there.'
	test ! -d $(FIG_DIR) && mkdir $(FIG_DIR) || echo '$(FIG_DIR)/ is already there.'

# get data
$(DATA)/%.tsv: $(SCRIPTS)/get_data.R
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# run the main part of the analysis
$(DATA)/rkip_wgcna.rda: $(SCRIPTS)/analysis.R
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# initiat files
initiat_files:
	for f in `ls $(FIG_SRC) | cut -d '.' -f1`; do echo test ! -e $(FIG_DIR)/$$f.png && touch $(FIG_DIR)/$$f.png; done;
	for f in `ls $(TAB_SRC) | cut -d '.' -f1`; do echo test ! -e $(TAB_DIR)/$$f.tex && touch $(TAB_DIR)/$$f.tex; done;

# generate tables
$(TAB_DIR)/%.tex: $(TAB_SRC)/%.R $(DATA)/rkip_wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate figures
$(FIG_DIR)/%.png: $(FIG_SRC)/%.R $(DATA)/rkip_wgcna.rda
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# generate workflow
#$(FIG_DIR)/workflow.png: $(FIG_SRC)/workflow.dot
#	dot -Tpng $(FIG_SRC)/workflow.dot -o $(FIG_DIR)/workflow.png

# clean up
.PHONY: clean
clean:
	rm -f *.pdf
	rm -f *.RData
