from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField
from wtforms.validators import DataRequired, Length

class DNAform(FlaskForm):
    downSeq=StringField("Enter downstream sequence here:", validators=[DataRequired(), Length(min=25)])
    upSeq = StringField("Enter upstream sequence here:", validators=[DataRequired(), Length(min=25)])
    mutation = StringField("Enter mutation:", validators=[DataRequired(), Length(min=1, max=1)])
    WT = StringField("Enter wildtype:", validators=[DataRequired(), Length(min=1, max=1)])
    submit=SubmitField('Enter')

class contactForm(FlaskForm):
    message = StringField("Leave us a message here:", validators=[DataRequired()])
    submit = SubmitField('Enter')