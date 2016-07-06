from app import app
from config import ENABLE_CACHE, SAMPLE_FASTA
from view import status_view, gallery_view
from view import prj_info_view, help_view
from view import index_view
from tool import render


@app.after_request
def add_header(response):
    """
    Add headers to both force latest IE rendering engine or Chrome Frame,
    and also to cache the rendered page for 10 minutes.
    """
    response.headers['X-UA-Compatible'] = 'IE=Edge,chrome=1'
    if ENABLE_CACHE == False:
        response.headers['Cache-Control'] = 'public, max-age=0'
    return response


@app.route('/about')
def about():
    return render('about.html', title='About')


@app.route('/sample_fasta')
def sample_fasta():
    with open(SAMPLE_FASTA) as f:
        sample = f.read()
        return sample
    

@app.errorhandler(404)
def page_not_found(e):
    return render('404.html'), 404

@app.errorhandler(500)
def page_not_found(e):
    return render('500.html'), 500