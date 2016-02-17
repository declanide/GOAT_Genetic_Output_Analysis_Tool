# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations


class Migration(migrations.Migration):

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='affy_s_rs',
            fields=[
                ('id_s', models.AutoField(serialize=False, primary_key=True)),
                ('id_rs', models.CharField(max_length=20, null=True, blank=True)),
                ('affy_5', models.CharField(max_length=45, null=True, blank=True)),
                ('affy_6', models.CharField(max_length=45, null=True, blank=True)),
            ],
            options={
                'db_table': 'affy_s_rs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='assoc',
            fields=[
                ('idassoc', models.AutoField(serialize=False, primary_key=True)),
                ('rs_id_assoc', models.CharField(max_length=45, null=True, blank=True)),
                ('chromosome', models.IntegerField(null=True, blank=True)),
                ('pos', models.BigIntegerField(null=True, blank=True)),
                ('allele_A', models.CharField(max_length=1, null=True, blank=True)),
                ('allele_B', models.CharField(max_length=1, null=True, blank=True)),
                ('av_max_post_call', models.FloatField(null=True, blank=True)),
                ('info_imp', models.FloatField(null=True, blank=True)),
                ('cohort_AA', models.FloatField(null=True, blank=True)),
                ('cohort_AB', models.FloatField(null=True, blank=True)),
                ('cohort_BB', models.FloatField(null=True, blank=True)),
                ('cohort_NULL', models.FloatField(null=True, blank=True)),
                ('cases_AA', models.FloatField(null=True, blank=True)),
                ('cases_AB', models.FloatField(null=True, blank=True)),
                ('cases_BB', models.FloatField(null=True, blank=True)),
                ('controls_AA', models.FloatField(null=True, blank=True)),
                ('controls_AB', models.FloatField(null=True, blank=True)),
                ('controls_BB', models.FloatField(null=True, blank=True)),
                ('allele_A_freq', models.FloatField(null=True, blank=True)),
                ('allele_B_freq', models.FloatField(null=True, blank=True)),
                ('maf', models.FloatField(null=True, blank=True)),
                ('cases_maf', models.FloatField(null=True, blank=True)),
                ('controls_maf', models.FloatField(null=True, blank=True)),
                ('miss_dat_prop', models.FloatField(null=True, blank=True)),
                ('cohort_hwe', models.FloatField(null=True, blank=True)),
                ('cases_hwe', models.FloatField(null=True, blank=True)),
                ('controls_hwe', models.FloatField(null=True, blank=True)),
                ('het_OR', models.FloatField(null=True, blank=True)),
                ('het_OR_lower', models.FloatField(null=True, blank=True)),
                ('het_OR_upper', models.FloatField(null=True, blank=True)),
                ('all_OR', models.FloatField(null=True, blank=True)),
                ('all_OR_lower', models.FloatField(null=True, blank=True)),
                ('all_OR_upper', models.FloatField(null=True, blank=True)),
                ('pvalue_assoc', models.FloatField(null=True, blank=True)),
                ('info_assoc', models.FloatField(null=True, blank=True)),
                ('beta_assoc', models.FloatField(null=True, blank=True)),
                ('se_assoc', models.FloatField(null=True, blank=True)),
            ],
            options={
                'db_table': 'assoc',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='covariate',
            fields=[
                ('idcovariate', models.AutoField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'covariate',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='dataset',
            fields=[
                ('iddataset', models.AutoField(serialize=False, primary_key=True)),
                ('name', models.CharField(max_length=45, null=True, blank=True)),
                ('description_dataset', models.CharField(max_length=45, null=True, blank=True)),
                ('cohort', models.CharField(max_length=45, null=True, blank=True)),
                ('cohort_description', models.CharField(max_length=45, null=True, blank=True)),
                ('build_id', models.IntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'dataset',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='experiment',
            fields=[
                ('idexperiment', models.AutoField(serialize=False, primary_key=True)),
                ('date_exp', models.DateTimeField(null=True, blank=True)),
                ('version_script', models.FloatField(null=True, blank=True)),
                ('verstion_algo', models.FloatField(null=True, blank=True)),
                ('post_process', models.TextField(null=True, blank=True)),
                ('note', models.TextField(null=True, blank=True)),
                ('covariates', models.CharField(max_length=45, null=True, blank=True)),
                ('experiment_approved', models.IntegerField()),
            ],
            options={
                'db_table': 'experiment',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='genes',
            fields=[
                ('idgenes', models.AutoField(serialize=False, primary_key=True)),
                ('symbol', models.CharField(max_length=45, null=True, blank=True)),
                ('approved', models.CharField(max_length=45, null=True, blank=True)),
                ('type_field', models.CharField(max_length=45, null=True, db_column='type', blank=True)),
                ('longname', models.CharField(max_length=150, null=True, blank=True)),
                ('note', models.TextField(null=True, blank=True)),
            ],
            options={
                'db_table': 'genes',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='hg_rs_history',
            fields=[
                ('idhg_rs_history', models.AutoField(serialize=False, primary_key=True)),
                ('rsHigh', models.CharField(max_length=20, null=True, blank=True)),
                ('rsLow', models.CharField(max_length=45, null=True, blank=True)),
                ('build_id', models.IntegerField(null=True, blank=True)),
                ('orien', models.BooleanField()),
                ('create_time', models.DateTimeField(null=True, blank=True)),
                ('last_updated_time', models.DateTimeField(null=True, blank=True)),
                ('rsCurrent', models.CharField(max_length=20, null=True, blank=True)),
                ('orien2Current', models.BooleanField()),
            ],
            options={
                'db_table': 'hg_rs_history',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='mar_gene',
            fields=[
                ('idmarq_gene', models.AutoField(serialize=False, primary_key=True)),
                ('nom', models.CharField(max_length=20, null=True, blank=True)),
                ('build_id', models.IntegerField()),
                ('symbol', models.CharField(max_length=50, null=True, blank=True)),
                ('pos', models.CharField(max_length=1, null=True, blank=True)),
                ('dist', models.IntegerField()),
                ('strand', models.CharField(max_length=1, null=True, blank=True)),
                ('anot', models.CharField(max_length=50, null=True, blank=True)),
            ],
            options={
                'db_table': 'covariate',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='marqueurs',
            fields=[
                ('idmarqueurs', models.AutoField(serialize=False, primary_key=True)),
                ('nom', models.CharField(unique=True, max_length=20, primary_key=True)),
                ('sorte', models.SmallIntegerField(null=True, blank=True)),
                ('chromosome', models.SmallIntegerField(null=True, blank=True)),
                ('position', models.BigIntegerField(null=True, blank=True)),
                ('build_Id', models.IntegerField(primary_key=True)),
                ('remapped_from_hg18', models.CharField(max_length=20, unique=True, null=True, blank=True)),
                ('refNCBI', models.CharField(max_length=5, null=True, blank=True)),
                ('observed', models.CharField(max_length=5, null=True, blank=True)),
                ('classe', models.CharField(max_length=15, null=True, blank=True)),
                ('func', models.CharField(max_length=15, null=True, blank=True)),
                ('frame', models.CharField(max_length=15, null=True, blank=True)),
                ('codons', models.CharField(max_length=10, null=True, blank=True)),
                ('peptides', models.CharField(max_length=5, null=True, blank=True)),
                ('gene', models.CharField(max_length=15, null=True, blank=True)),
                ('gene_strand', models.CharField(max_length=15, null=True, blank=True)),
                ('start_gen', models.BigIntegerField(null=True, blank=True)),
                ('end_gen', models.BigIntegerField(null=True, blank=True)),
                ('gene_before', models.CharField(max_length=15, null=True, blank=True)),
                ('gene_before_strand', models.CharField(max_length=15, null=True, blank=True)),
                ('dist_gen_before', models.IntegerField(null=True, blank=True)),
                ('start_gen_before', models.BigIntegerField(null=True, blank=True)),
                ('end_gen_before', models.BigIntegerField(null=True, blank=True)),
                ('gene_after', models.CharField(max_length=15, null=True, blank=True)),
                ('gene_after_strand', models.CharField(max_length=15, null=True, blank=True)),
                ('dist_gen_after', models.IntegerField(null=True, blank=True)),
                ('start_gen_after', models.BigIntegerField(null=True, blank=True)),
                ('end_gen_after', models.BigIntegerField(null=True, blank=True)),
            ],
            options={
                'db_table': 'marqueurs',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='person',
            fields=[
                ('idperson', models.AutoField(serialize=False, primary_key=True)),
                ('sex', models.CharField(max_length=15, null=True, blank=True)),
                ('age', models.IntegerField(null=True, blank=True)),
                ('dob', models.DateTimeField(null=True, blank=True)),
                ('ethnic_code', models.IntegerField(null=True, blank=True)),
                ('country_name', models.CharField(max_length=45, null=True, blank=True)),
                ('region_name', models.CharField(max_length=45, null=True, blank=True)),
                ('centre_name', models.CharField(max_length=200, null=True, blank=True)),
                ('sample', models.IntegerField(null=True, blank=True)),
                ('comment', models.TextField(null=True, blank=True)),
            ],
            options={
                'db_table': 'person',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='person_dataset',
            fields=[
                ('idperson_dataset', models.IntegerField(serialize=False, primary_key=True)),
            ],
            options={
                'db_table': 'person_dataset',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='phen_id_val',
            fields=[
                ('idphen_id_val', models.AutoField(serialize=False, primary_key=True)),
                ('value', models.FloatField(null=True, blank=True)),
            ],
            options={
                'db_table': 'phen_id_val',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='phenotypes',
            fields=[
                ('idphenotypes', models.AutoField(serialize=False, primary_key=True)),
                ('nom', models.CharField(max_length=45, null=True, blank=True)),
                ('nom_view', models.CharField(max_length=45, null=True, blank=True)),
                ('description', models.TextField(null=True, blank=True)),
                ('type_field', models.CharField(max_length=45, null=True, db_column='type', blank=True)),
                ('querry', models.TextField(null=True, blank=True)),
                ('covariable_choosing_table', models.CharField(max_length=200, null=True, blank=True)),
                ('characteristic_table', models.CharField(max_length=200, null=True, blank=True)),
            ],
            options={
                'db_table': 'phenotypes',
                'managed': False,
            },
        ),
        migrations.CreateModel(
            name='users',
            fields=[
                ('idusers', models.AutoField(serialize=False, primary_key=True)),
                ('nom', models.CharField(max_length=45, null=True, blank=True)),
                ('note', models.CharField(max_length=45, null=True, blank=True)),
                ('password', models.CharField(max_length=45, null=True, blank=True)),
            ],
            options={
                'db_table': 'users',
                'managed': False,
            },
        ),
    ]
