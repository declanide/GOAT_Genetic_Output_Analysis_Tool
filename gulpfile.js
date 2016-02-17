var gulp = require('gulp');
var shell = require('gulp-shell');
var changed = require('gulp-changed');
var plumber = require('gulp-plumber');
var sass = require ('gulp-sass');

var OPTIONS = {
    COLLECSTATIC:{
        watch:'static/**/*.*'
    },
    SASS:{
        src:'static/css/scss/app.scss',
        dest:'static/css'
    }
}


gulp.task('sass',function(){
  return gulp.src('static/css/scss/app.scss')
    .pipe(sass().on('error',sass.logError))
    .pipe(gulp.dest('static/css'))
});

gulp.task('collectstatic',['sass'], shell.task([
  'python manage.py collectstatic --noinput --ignore "*.scss"'
]));

gulp.task('watch', function () {
    gulp.watch('static/**/*.*');
    gulp.watch('biomarqueurs/templates/**/*.*');

});

gulp.task('default',['collectstatic'],function(){

})