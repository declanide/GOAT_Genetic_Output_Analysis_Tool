from django.contrib import admin
from models import *

class assocAdmin(admin.ModelAdmin):
    list_display = ('idassoc','experiment', 'rs_id_assoc', 'chromosome' ,'pos','allele_A', 'allele_B', 'av_max_post_call','info_imp','cohort_AA','cohort_AB','cohort_BB','cohort_NULL','cases_AA', 'cases_AB', 'cases_BB', 'controls_AA', 'controls_AB', 'controls_BB', 'allele_A_freq', 'allele_B_freq', 'maf', 'cases_maf', 'controls_maf', 'miss_dat_prop', 'cohort_hwe', 'cases_hwe', 'controls_hwe', 'het_OR', 'het_OR_lower', 'het_OR_upper', 'all_OR','all_OR_lower','all_OR_upper', 'pvalue_assoc', 'info_assoc','beta_assoc', 'se_assoc')
    search_fields = ['idassoc','experiment', 'rs_id_assoc', 'chromosome' ,'pos','allele_A', 'allele_B', 'av_max_post_call','info_imp','cohort_AA','cohort_AB','cohort_BB','cohort_NULL','cases_AA', 'cases_AB', 'cases_BB', 'controls_AA', 'controls_AB', 'controls_BB', 'allele_A_freq', 'allele_B_freq', 'maf', 'cases_maf', 'controls_maf', 'miss_dat_prop', 'cohort_hwe', 'cases_hwe', 'controls_hwe', 'het_OR', 'het_OR_lower', 'het_OR_upper', 'all_OR','all_OR_lower','all_OR_upper', 'pvalue_assoc', 'info_assoc','beta_assoc', 'se_assoc']

class marqueursAdmin(admin.ModelAdmin):
    list_display = ('idmarqueurs', 'nom', 'sorte','chromosome','position', 'idgenes','build_Id', 'remapped_from_hg18', 'refNCBI', 'observed', 'classe', 'func', 'frame', 'codons', 'peptides', 'gene', 'gene_strand', 'start_gen', 'end_gen', 'gene_before', 'gene_before_strand','dist_gen_before', 'start_gen_before', 'end_gen_before', 'gene_after', 'gene_after_strand', 'dist_gen_after', 'start_gen_after', 'end_gen_after')
    search_fields = ['idmarqueurs', 'nom', 'sorte','chromosome','position', 'idgenes','build_Id', 'remapped_from_hg18', 'refNCBI', 'observed', 'classe', 'func', 'frame', 'codons', 'peptides', 'gene', 'gene_strand', 'start_gen','end_gen', 'gene_before', 'gene_before_strand','dist_gen_before', 'start_gen_before', 'end_gen_before', 'gene_after', 'gene_after_strand', 'dist_gen_after', 'start_gen_after', 'end_gen_after']

class hg_rs_historyAdmin(admin.ModelAdmin):
    list_display = ('idhg_rs_history', 'rsHigh', 'rsLow','build_id','orien','create_time','last_updated_time','rsCurrent', 'orien2Current')
    search_fields = ['idhg_rs_history ', 'rsHigh', 'rsLow','build_id','orien','create_time','last_updated_time','rsCurrent', 'orien2Current']

class affy_s_rsAdmin(admin.ModelAdmin):
    list_display = ('id_s','id_rs', 'affy_5', 'affy_6')
    search_fields = ['id_s','id_rs', 'affy_5', 'affy_6']

class genesAdmin(admin.ModelAdmin):
    list_display = ('idgenes', 'symbol','approved','type_field','longname','note')
    search_fields= ['idgenes', 'symbol','approved','type_field','longname','note']

class experimentAdmin(admin.ModelAdmin):
    list_display = ('idexperiment', 'idusers','date_exp', 'version_script','verstion_algo','post_process','note','dataset','phenotype','covariates','experiment_approved')
    search_fields = ['idexperiment', 'idusers','date_exp', 'version_script','verstion_algo','post_process','note','dataset','phenotype','covariates','experiment_approved']

class usersAdmin(admin.ModelAdmin):
    list_display = ('idusers','nom','note','password')
    search_fields = ['idusers','nom','note','password']

class phenotypesAdmin(admin.ModelAdmin):
    list_display = ('idphenotypes','nom','nom_view','description','type_field','querry','covariable_choosing_table', 'characteristic_table')
    search_fields = ['idphenotypes','nom','nom_view','description','type_field','querry','covariable_choosing_table', 'characteristic_table']

class phen_id_valAdmin(admin.ModelAdmin):
    list_display = ('idphen_id_val','idphenotypes','idperson','value')
    search_fields = ['idphen_id_val','idphenotypes','idperson','value']

class datasetAdmin(admin.ModelAdmin):
    list_display = ('iddataset','name','description_dataset', 'cohort', 'cohort_description','build_id')
    search_fields = ['iddataset','name','description_dataset', 'cohort', 'cohort_description','build_id']

class person_datasetAdmin(admin.ModelAdmin):
    list_display = ('idperson_dataset','idperson', 'iddataset')
    search_fields = ['idperson_dataset','idperson', 'iddataset']

class personAdmin(admin.ModelAdmin):
    list_display = ('idperson', 'sex','age', 'dob','ethnic_code','country_name','region_name', 'centre_name', 'sample', 'comment')
    search_fields = ['idperson', 'sex','age', 'dob','ethnic_code','country_name','region_name', 'centre_name', 'sample', 'comment']

class covariateAdmin(admin.ModelAdmin):
    list_display = ('idcovariate','idexperiment','idphenotypes')
    search_fields = ['idcovariate','idexperiment','idphenotypes']

class marq_geneAdmin(admin.ModelAdmin):
    list_display = ('idmarq_gene','idmarqueurs','nom','build_id','idgenes','symbol','pos','dist','strand','anot')
    search_fields = ['idmarq_gene','idmarqueurs','nom','build_id','idgenes','symbol','pos','dist','strand','anot']

admin.site.register(mar_gene, marq_geneAdmin)
admin.site.register(covariate,covariateAdmin)
admin.site.register(person,personAdmin)
admin.site.register(person_dataset,person_datasetAdmin)
admin.site.register(dataset,datasetAdmin)
admin.site.register(phen_id_val,phen_id_valAdmin)
admin.site.register(users,usersAdmin)
admin.site.register(phenotypes,phenotypesAdmin)
admin.site.register(experiment, experimentAdmin)
admin.site.register(genes, genesAdmin)
admin.site.register(affy_s_rs,affy_s_rsAdmin)
admin.site.register(assoc,assocAdmin)
admin.site.register(hg_rs_history,hg_rs_historyAdmin)
admin.site.register(marqueurs,marqueursAdmin)



