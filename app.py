from flask import Flask, render_template, url_for, flash, redirect,request,jsonify
from forms import DNAform, contactForm,crisPAMform, uploadForm
from Scripts.BEsingle import MainBE
from Scripts.seqAnalysis import process_sequence
import os
import secrets

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

def save_file(input_field):
    random_hex=secrets.token_hex(8)
    f_name, f_ext=os.path.splitext(input_field.filename)
    filename=random_hex+f_ext
    filepath=os.path.join(app.root_path, "user_files",filename)
    input_field.save(filepath)


@app.route('/upload',methods=['GET', 'POST'])
def upload():
    form=uploadForm()
    result=''
    if form.validate_on_submit():
        save_file(form.inputFile.data)
        result="123"
    elif request.method=='GET':
        result='1'

    return render_template('upload.html',form=form,result=result)

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
