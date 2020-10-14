from os import path
import psutil

### environment variables ###
ENVIR = {
'PATH' : '/usr/bin/:\
/home/xuyihui/RGAugury/ncbi-blast-2.10.1+-x64-linux/ncbi-blast-2.10.1+/bin:\
/home/xuyihui/RGAugury/hmmer3.1b2:\
/home/xuyihui/RGAugury/PfamScan:\
/home/xuyihui/database/interproscan-5.46-81.0-64-bit/interproscan-5.46-81.0:\
/home/xuyihui/RGAugury/rgaugury:\
/home/xuyihui/RGAugury/rgaugury/coils:\
/home/xuyihui/RGAugury/cvit.1.2.1:\
/home/xuyihui/RGAugury/phobius101_linux/tmp/tmprnL5On/phobius',

 'JAVA_HOME':'/usr/lib/jvm/java-11-openjdk-amd64',
 'PERL5LIB':'/home/xuyihui/RGAugury/PfamScan',
 'COILSDIR':'/home/xuyihui/RGAugury/rgaugury/coils',
 'PFAMDB':'/home/xuyihui/database/pfamdb'
}

### production environment ########
WEB_UI_LOG = False
ENABLE_CACHE = True      ### ! important, set enable_cache as True in the deployment enviroment.
DEBUG = False   ### ! important
PIPELINE_HOME='/home/xuyihui/RGAugury/rgaugury'
# any bitbucket username and password will work
USER_NAME = 'xuyihui1999@icloud.com'
PASSWORD = 'Xu65183843'

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
