from app import db

class User(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    username = db.Column(db.String(64),unique=True)
    password = db.Column(db.String(64))
    email = db.Column(db.String(120),unique=True)

class GeneAmount(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    nbs = db.Column(db.Integer)
    cnl = db.Column(db.Integer)
    tnl = db.Column(db.Integer)
    cn = db.Column(db.Integer)
    tn = db.Column(db.Integer)
    nl = db.Column(db.Integer)
    tx = db.Column(db.Integer)
    other = db.Column(db.Integer)
    rlp = db.Column(db.Integer)
    rlk = db.Column(db.Integer)
    tmcc = db.Column(db.Integer)
    total = db.Column(db.Integer)


class Project(db.Model):
    id = db.Column(db.String(64), primary_key = True)
    name = db.Column(db.String(50))
    status = db.Column(db.String(120))
    start_time = db.Column(db.String(64))
    elapsed_time = db.Column(db.String(64))
    gff3 = db.Column(db.Integer)
    e_value = db.Column(db.String(32))
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
    gene_amount = db.Column(db.Integer, db.ForeignKey('gene_amount.id'))
    pid = db.Column(db.Integer)
    pid_ctime = db.Column(db.Float)
    fingerprint = db.Column(db.Integer)
    step = db.Column(db.String(1000))
    input_amount = db.Column(db.Integer)
    database = db.Column(db.String(50))

class Gene(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    prj_id = db.Column(db.Integer, db.ForeignKey('project.id'))
    name = db.Column(db.String(64))
    length = db.Column(db.Integer)
    type = db.Column(db.String(32))
    image = db.Column(db.String(120))

    def to_list(self):
        types = self.type.split('>')
        return [self.name, self.length, types[-1], self.image]

class GeneSection(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    span = db.Column(db.String(20))
    gene_id =db.Column(db.Integer, db.ForeignKey('gene.id'))
    category = db.Column(db.String(20))
    def toList(self):
        return [self.span, self.category]

class GeneImage(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(20))
    category = db.Column(db.String(20))

class Gff(db.Model):
    seqname = db.Column(db.String(20), primary_key=True)
    source = db.Column(db.String(20), primary_key=True)
    feature = db.Column(db.String(20), primary_key=True)
    start = db.Column(db.String(20), primary_key=True)
    end = db.Column(db.String(20), primary_key=True)
    score = db.Column(db.String(20), primary_key=True)
    strand = db.Column(db.String(20), primary_key=True)
    frame = db.Column(db.String(20), primary_key=True)
    attribute  = db.Column(db.String(20), primary_key=True)