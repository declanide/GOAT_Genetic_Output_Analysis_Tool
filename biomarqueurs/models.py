# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# Beatriz Kanzki
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# -----------------------------------------------------------Modeles de la base de donnees hits-of-opti-tera cree pour logiciel "GOAT" avec Django--------------------------------------------------------------#
# -------------------------------------------------------------------------------------------13/05/2015---------------------------------------------------------------------------------------------------------#
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
from __future__ import unicode_literals

from django.db import models

class ResultsManager(models.Manager):
    def get_assoc(self):
        print("test2")
        return self


class assoc(models.Model):
    idassoc = models.AutoField(primary_key=True)
    experiment = models.ForeignKey('biomarqueurs.experiment', db_column='experiment', blank=True, null=True)
    rs_id_assoc = models.CharField(max_length=45, blank=True, null=True)
    chromosome = models.IntegerField(blank=True, null=True)
    pos = models.BigIntegerField(blank=True, null=True)
    allele_A = models.CharField(max_length=1, blank=True, null=True)
    allele_B = models.CharField(max_length=1, blank=True, null=True)
    av_max_post_call = models.FloatField(blank=True, null=True)
    info_imp = models.FloatField(blank=True, null=True)
    cohort_AA = models.FloatField(blank=True, null=True)
    cohort_AB = models.FloatField(blank=True, null=True)
    cohort_BB = models.FloatField(blank=True, null=True)
    cohort_NULL = models.FloatField(blank=True, null=True)
    cases_AA = models.FloatField(blank=True, null=True)
    cases_AB = models.FloatField(blank=True, null=True)
    cases_BB = models.FloatField(blank=True, null=True)
    controls_AA = models.FloatField(blank=True, null=True)
    controls_AB = models.FloatField(blank=True, null=True)
    controls_BB = models.FloatField(blank=True, null=True)
    allele_A_freq = models.FloatField(blank=True, null=True)
    allele_B_freq = models.FloatField(blank=True, null=True)
    maf = models.FloatField(blank=True, null=True)
    cases_maf = models.FloatField(blank=True, null=True)
    controls_maf = models.FloatField(blank=True, null=True)
    miss_dat_prop = models.FloatField(blank=True, null=True)
    cohort_hwe = models.FloatField(blank=True, null=True)
    cases_hwe = models.FloatField(blank=True, null=True)
    controls_hwe = models.FloatField(blank=True, null=True)
    het_OR = models.FloatField(blank=True, null=True)
    het_OR_lower = models.FloatField(blank=True, null=True)
    het_OR_upper = models.FloatField(blank=True, null=True)
    all_OR = models.FloatField(blank=True, null=True)
    all_OR_lower = models.FloatField(blank=True, null=True)
    all_OR_upper = models.FloatField(blank=True, null=True)
    pvalue_assoc = models.FloatField(blank=True, null=True)
    info_assoc = models.FloatField(blank=True, null=True)
    beta_assoc = models.FloatField(blank=True, null=True)
    se_assoc = models.FloatField(blank=True, null=True)
    objects = ResultsManager()
    class Meta:
        managed = True
        db_table = 'assoc'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class marqueurs(models.Model):
    idmarqueurs = models.AutoField(primary_key=True)
    nom = models.CharField(max_length= 20, unique = True, primary_key=True)
    sorte = models.SmallIntegerField(blank=True, null=True)
    chromosome = models.SmallIntegerField(blank=True, null=True)
    position = models.BigIntegerField(blank=True, null=True)
    idgenes = models.ForeignKey('biomarqueurs.genes', db_column='idgenes', blank=True, null=True)
    build_Id = models.IntegerField(primary_key=True)
    remapped_from_hg18 = models.CharField(unique=True, max_length=20, blank=True, null=True)
    refNCBI = models.CharField(max_length=5, blank=True, null=True)
    observed = models.CharField(max_length=5, blank=True, null=True)
    classe = models.CharField(max_length=15, blank=True, null=True)
    func = models.CharField(max_length=15, blank=True, null=True)
    frame = models.CharField(max_length=15, blank=True, null=True)
    codons = models.CharField(max_length=10, blank=True, null=True)
    peptides = models.CharField(max_length=5, blank=True, null=True)
    gene = models.CharField(max_length=15, blank=True, null=True)
    gene_strand = models.CharField(max_length=15, blank=True, null=True)
    start_gen = models.BigIntegerField(blank=True, null=True)
    end_gen = models.BigIntegerField(blank=True, null=True)
    gene_before = models.CharField(max_length=15, blank=True, null=True)
    gene_before_strand = models.CharField(max_length=15, blank=True, null=True)
    dist_gen_before = models.IntegerField(blank=True, null=True)
    start_gen_before = models.BigIntegerField(blank=True, null=True)
    end_gen_before = models.BigIntegerField(blank=True, null=True)
    gene_after = models.CharField(max_length=15, blank=True, null=True)
    gene_after_strand = models.CharField(max_length=15, blank=True, null=True)
    dist_gen_after = models.IntegerField(blank=True, null=True)
    start_gen_after = models.BigIntegerField(blank=True, null=True)
    end_gen_after = models.BigIntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'marqueurs'
        unique_together = (('nom', 'remapped_from_hg18'),)


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class hg_rs_history(models.Model):
    idhg_rs_history = models.AutoField(primary_key=True)
    rsHigh = models.CharField(max_length=20, blank=True, null=True)
    rsLow = models.CharField(max_length=45, blank=True, null=True)
    build_id = models.IntegerField(blank=True, null=True)
    orien = models.BooleanField()
    create_time = models.DateTimeField(blank=True, null=True)
    last_updated_time = models.DateTimeField(blank=True, null=True)
    rsCurrent = models.CharField(max_length=20, blank=True, null=True)
    orien2Current = models.BooleanField()

    class Meta:
        managed = False
        db_table = 'hg_rs_history'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class affy_s_rs(models.Model):
    id_s = models.AutoField(primary_key=True)
    id_rs = models.CharField(max_length=20, blank=True, null=True)
    affy_5 = models.CharField(max_length=45, blank=True, null=True)
    affy_6 = models.CharField(max_length=45, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'affy_s_rs'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class genes(models.Model):
    idgenes = models.AutoField(primary_key=True)
    symbol = models.CharField(max_length=45, blank=True, null=True)
    approved = models.CharField(max_length=45, blank=True, null=True)
    type_field = models.CharField(db_column='type', max_length=45, blank=True, null=True)
    longname = models.CharField(max_length=150, blank=True, null=True)
    note = models.TextField(blank=True, null=True)

    class Meta :
        managed = False
        db_table = 'genes'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class experiment(models.Model):
    idexperiment = models.AutoField(primary_key=True)
    idusers = models.ForeignKey('biomarqueurs.users', db_column= 'idusers', blank=True, null=True)
    date_exp = models.DateTimeField(blank=True, null=True)
    version_script = models.FloatField( blank=True, null=True)
    verstion_algo = models.FloatField( blank=True, null=True)
    post_process = models.TextField(blank=True, null=True)
    note = models.TextField(blank=True, null=True)
    dataset=models.ForeignKey('biomarqueurs.dataset', db_column= 'iddataset', blank=True, null=True)
    phenotype=models.ForeignKey('biomarqueurs.phenotypes', db_column= 'idphenotypes', blank=True, null=True)
    covariates=models.CharField(max_length=45, blank=True, null=True)
    experiment_approved=models.IntegerField()

    class Meta:
        managed = False
        db_table = 'experiment'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class users(models.Model):
    idusers = models.AutoField(primary_key=True)
    nom = models.CharField(max_length=45, blank=True, null=True)
    note = models.CharField(max_length=45, blank=True, null=True)
    password = models.CharField(max_length=45, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'users'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class phenotypes(models.Model):
    idphenotypes = models.AutoField(primary_key=True)
    nom = models.CharField(max_length=45, blank=True, null=True)
    nom_view = models.CharField(max_length=45, blank=True, null=True)
    description = models.TextField(blank=True, null=True)
    type_field = models.CharField(db_column='type', max_length=45, blank=True, null=True)
    querry = models.TextField(blank=True, null=True)
    covariable_choosing_table = models.CharField(max_length =200, blank=True, null=True)
    characteristic_table = models.CharField(max_length=200, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'phenotypes'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class phen_id_val(models.Model):
    idphen_id_val = models.AutoField(primary_key=True)
    idphenotypes = models.ForeignKey('biomarqueurs.phenotypes', db_column= 'idphenotypes', blank=True, null=True)
    idperson = models.ForeignKey('biomarqueurs.person', db_column= 'idperson')
    value = models.FloatField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'phen_id_val'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class dataset(models.Model):
    iddataset = models.AutoField(primary_key=True)
    name = models.CharField(max_length=45, blank=True, null=True)
    description_dataset = models.CharField(max_length=45, blank=True, null=True)
    cohort = models.CharField(max_length=45, blank=True, null=True)
    cohort_description = models.CharField(max_length=45, blank=True, null=True)
    build_id=models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'dataset'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class person_dataset(models.Model):
    idperson_dataset = models.IntegerField(primary_key=True)
    idperson = models.ForeignKey('biomarqueurs.person', db_column='idperson', blank=True, null=True)
    iddataset = models.ForeignKey('biomarqueurs.dataset', db_column='iddataset', blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'person_dataset'


# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class person(models.Model):
    idperson = models.AutoField(primary_key=True)
    sex = models.CharField(max_length=15, blank=True, null=True)
    age = models.IntegerField(blank=True, null=True)
    dob = models.DateTimeField(blank=True, null=True)
    ethnic_code = models.IntegerField(blank=True, null=True)
    country_name = models.CharField(max_length=45, blank=True, null=True)
    region_name = models.CharField(max_length=45, blank=True, null=True)
    centre_name = models.CharField(max_length=200, blank=True, null=True)
    sample = models.IntegerField(blank=True, null=True)
    comment = models.TextField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'person'
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
class covariate(models.Model):
    idcovariate = models.AutoField(primary_key=True)
    idexperiment = models.ForeignKey('biomarqueurs.experiment', db_column='idexperiment', blank=True, null=True)
    idphenotypes = models.ForeignKey('biomarqueurs.phenotypes',db_column='idphenotypes', blank=True, null=True)

    class Meta:
        managed= False
        db_table = 'covariate'
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

class mar_gene(models.Model):
    idmarq_gene = models.AutoField(primary_key=True)
    idmarqueurs = models.ForeignKey('biomarqueurs.marqueurs', db_column='idmarqueurs', blank=True, null=True)
    nom=models.CharField(max_length=20, blank=True, null=True)
    build_id=models.IntegerField()
    idgenes = models.ForeignKey('biomarqueurs.genes',db_column='idgenes', blank=True, null=True)
    symbol=models.CharField(max_length=50, blank=True, null=True)
    pos=models.CharField(max_length=1, blank=True, null=True)
    dist=models.IntegerField()
    strand=models.CharField(max_length=1, blank=True, null=True)
    anot=models.CharField(max_length=50, blank=True, null=True)

    class Meta:
        managed= False
        db_table = 'covariate'
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
# END
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

