from app import app, db, models
from config import PRJ_HOME, META_FILE, MOTIF_PATH, GENE_IMAGE
from app.tool import render, getGeneId, getSequence, getGff
import sys
from flask import jsonify,send_from_directory, url_for

class Span:
    start=''
    ending=''
    img=''
    category=''
    def __init__(self, start, ending, img,category):
        self.start = start
        self.ending = ending
        self.image = img
        self.category = category

@app.route('/gallery/img/<prj_id>/img/<img>')
def gallery(prj_id,img):
  img_path = url_for('index')+'img/'+prj_id+'/img/'+img
  name = '.'.join(img.split('.')[:-1])
  gene = db.session.query(models.Gene).filter(models.Gene.prj_id == prj_id,models.Gene.name == name).first()
  gene_type = gene.type.split('>')[-1]
  sections = getSections(prj_id,name)
  spans = []
  for section in sections:
    geneImage = db.session.query(models.GeneImage).filter_by(category = section.category).first()
    l = section.span.split('|')
    start = l[0]
    ending = l[1]
    spans.append(Span(start,ending,url_for('index')+GENE_IMAGE +'/'+geneImage.name,section.category))
  
  return render('gallery.html',root = url_for('index'), title='GSV', prj_id=prj_id, img_path=img_path, gffs= getGff(prj_id,gene.name), \
                               name= gene.name, type =gene_type, length = gene.length, simpleName =gene.name.split('.')[0], \
                               sections =spans, sequence = getSequence(prj_id, gene.name))


@app.route('/section/<prj_id>/<geneName>')
def jsonifiedSection(prj_id,geneName):
    sections = getSections(prj_id,geneName)
    return jsonify(data=[section.toList() for section in sections])

@app.route('/'+GENE_IMAGE+'/<name>')
def getGeneImage(name):
    return send_from_directory(MOTIF_PATH,name)

def getSections(prj_id,geneName):
    gene_id =  getGeneId(prj_id,geneName)
    sections = db.session.query(models.GeneSection).filter_by(gene_id = gene_id).all()
    return sections