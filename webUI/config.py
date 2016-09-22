from os import path
import psutil

### environment variables ###
ENVIR = {
'PATH' : '/home/quanx/app/jdk1.8.0_91/bin:/home/quanx/app/ncbi-blast-2.3.0+/bin:\
/home/quanx/app/hmmer-3.1b2-linux-intel-x86_64/binaries:/home/quanx/app/PfamScan:\
/home/quanx/app/interproscan-5.17-57.0:/home/quanx/repo/rga:\
/home/quanx/repo/rga/coils:/home/quanx/app/cvit.1.2.1:/home/quanx/app/phobius1.01:\
/home/quanx/app/ncbi-blast-2.3.0+/bin',

'JAVA_HOME':'/home/quanx/app/jdk1.8.0_91',
'PERL5LIB':'/home/quanx/app/PfamScan',
'COILSDIR':'/home/quanx/repo/rga/coils',
'PFAMDB':'/home/quanx/pfamdb'
}

### production environment ########
WEB_UI_LOG = False
ENABLE_CACHE = False       ### ! important, set enable_cache as True in the deployment enviroment.
DEBUG = True   ### ! important
PIPELINE_HOME='/home/quanx/repo/rga'
# any bitbucket username and password will work
USER_NAME = ''
PASSWORD = ''

### application ###
PORT = 7000
CPU_TOGGLE = 0 # 1 means that users can choose cpu amount used to run the pipeline
APP_HOME = 'app' # Do not change this value
BASE_PATH = path.abspath(path.dirname(__file__))
TEMPLATE = 'templates'
### database ###
DB_DIR = path.join(BASE_PATH, APP_HOME, 'database')
DB_PATH = path.join(DB_DIR, 'app.db')
SQLALCHEMY_DATABASE_URI = 'sqlite:///' + DB_PATH
SQLALCHEMY_MIGRATE_REPO = path.join(DB_DIR, 'db_repository')

### website ###
ADMIN_PASSWORD = '3355'
WTF_CSRF_ENABLED = True
SECRET_KEY = 'you-will-never-guess'
PROJECTS ='/projects'
PRJ_IMG_PATH = '/img'
PRJ_HOME = BASE_PATH + PROJECTS
MOTIF_PATH = PIPELINE_HOME + '/motif'
PERL='perl'
RGAUGURY_PL='RGAugury.pl'

### project ###
TOTAL_STEPS = 15
FINAL_STATUS = 'complete' ## If you change this string, you also need to change code in status.html
CANCELED ='canceled' ## If you change this string, you also need to change code in status.html
PENDING ='pending'
FASTA_EXTENSION = '.RGA.candidates.fasta'
GENE_FILE = '.RGA.info.txt'
SUMMARY_FILE = '.RGA.summaries.txt'
META_FILE ='.RGA.gene.meta.txt'
START_PIPELINE = BASE_PATH + '/startPipeline.py'
LOG_FILE = '.status.log'
GFF_FILE = '.RGA.gff'
WEB_LOG = 'webui.log'
GENE_IMAGE='gene_image'
SAMPLE_GFF = 'sample.gff'
### others ###
DATE_FORMAT = '%Y/%m/%d %H:%M:%S'
SAMPLE_FASTA = BASE_PATH + '/sample.fas'
TOTAL_CPU = psutil.cpu_count()
JOBS = 1 # max number of jobs that run concurrently

### wiki help ###
WIKI = 'https://bitbucket.org/yaanlpc/rgaugury/wiki/Web%20UI%20Help'
