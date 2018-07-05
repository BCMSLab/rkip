R_OPT=--vanilla

SCRIPTS=scripts
FIG_SRC=scripts/figures
TAB_SRC=scripts/tables

MANUSCRIPT=manuscript
FIG_DIR=manuscript/figures
TAB_DIR=manuscript/tables

DATA=data
LOG=log

all: dir_tree \
	$(LOG)/get_data.R.Rout \
	$(DATA)/rkip_wgcna.rda \
	$(FIG_DIR)/colocalization.png \
	$(FIG_DIR)/colocalization_supp.png \
	$(FIG_DIR)/common_regulation.png \
	$(FIG_DIR)/gene_membership.png \
	$(FIG_DIR)/interactions_profile.png \
	$(FIG_DIR)/module_cor.png \
	$(FIG_DIR)/network.png \
	$(FIG_DIR)/pebp1_level.png \
	$(FIG_DIR)/preservation_ranks.png \
	$(FIG_DIR)/preservation_zsummary.png \
	$(FIG_DIR)/sft_power.png \
	$(TAB_DIR)/datasets.tex \
	$(TAB_DIR)/go_genes.tex \
	$(TAB_DIR)/module_members.tex \
	$(TAB_DIR)/peb_interactions.tex \
	$(TAB_DIR)/previously_reported.tex \
	clean

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
$(LOG)/get_data.R.Rout: $(SCRIPTS)/get_data.R
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

# run the main part of the analysis
$(DATA)/rkip_wgcna.rda: $(SCRIPTS)/analysis.R
	R CMD BATCH $(R_OPT) $< $(LOG)/$(<F).Rout

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
