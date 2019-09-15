from flask import Flask, render_template, url_for, flash, redirect,request,jsonify
from forms import DNAform, contactForm,crisPAMform
from Scripts.BEsingle import MainBE
from Scripts.seqAnalysis import process_sequence
import json

app = Flask(__name__)
app.config['SECRET_KEY']="3af987093f78a156771d3e00335bd60c"
posts = [
    {
        'author': 'Shiri Almog',
        'title': 'Blog Post 1',
        'content': 'First post content',
        'date_posted': 'Sept 10, 2019'
    },
    {
        'author': '123',
        'title': 'Blog Post 2',
        'content': 'Second post content',
        'date_posted': 'April 21, 2018'
    }
]

@app.route('/')
@app.route('/home')
def hello_world():
    return render_template('home.html', posts=posts)
@app.route('/blog')
def hello_world2():
    return render_template('blog.html', posts=posts)
@app.route('/about')
def about_world():
    return render_template('about.html')
@app.route('/contact')
def contact():
    form=contactForm()
    return render_template('contact.html',title="contact", form=form)


@app.route('/upload',methods=['GET', 'POST'])
def upload():
    return render_template('upload.html')

@app.route('/crisPAM', methods=['GET','POST'])
def crisPAM():
    form=crisPAMform()
    result=''
    if form.validate_on_submit():
        result=process_sequence(form.referencePAM.data,form.variantPAM.data)
    return render_template('crisPAM.html',title='crisPAM',form=form,result=result)


@app.route('/input', methods=['GET', 'POST'])
def input():
    form=DNAform()
    matches = ''
    clean = ''
    quiet = ''
    readingFrame = 2
    showReadingFrame = False
    if form.validate_on_submit():
        matches,clean,quiet = MainBE(form.upSeq.data,form.downSeq.data,form.mutation.data,form.WT.data)
    return render_template('input.html', title="input", form=form,
                           result=(matches,clean,quiet), readingFrame=readingFrame, showReadingFrame=showReadingFrame)

if __name__ == '__main__':
    app.run(debug=True)
