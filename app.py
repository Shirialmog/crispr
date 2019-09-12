from flask import Flask, render_template, url_for, flash, redirect,request,jsonify
from forms import DNAform, contactForm
from Scripts.process import myadd
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
@app.route('/process', methods=['post'])
def process():
    name=request.form['message']
    return myadd(name)

@app.route('/input', methods=['GET', 'POST'])
def input():
    form=DNAform()
    myresult = ''
    if form.validate_on_submit():
        myresult = myadd(form.downSeq.data)
        #return redirect(url_for('hello_world'))
    return render_template('input.html', title="inputzzz", form=form, result=myresult)

if __name__ == '__main__':
    app.run(debug=True)
