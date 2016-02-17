"""goat URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.8/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Add an import:  from blog import urls as blog_urls
    2. Add a URL to urlpatterns:  url(r'^blog/', include(blog_urls))
"""
from django.conf.urls import patterns, include, url
from django.contrib import admin

from biomarqueurs.views import *
from django.conf.urls import include, url
from django.conf import settings
from GOAT_Genetic_Ouput_Analysis_Tool import settings
from django.conf.urls.static import static
from django.contrib import admin
from biomarqueurs.views import index
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

admin.autodiscover()

urlpatterns = [
    url(r'^admin/', include(admin.site.urls)),
    url(r'^login/','biomarqueurs.views.login',name="login"),
    url(r'^logout/','biomarqueurs.views.logoutUser',name="logoutUser"),
    url(r'^SGene/','biomarqueurs.views.SGene',name="search"),
    url(r'^','biomarqueurs.views.index',name="home"),

]+static(settings.STATIC_URL,document_root=settings.STATIC_ROOT)+ static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) #done only in development not in production

if settings.DEBUG:              #only for development
    urlpatterns+=staticfiles_urlpatterns()
    urlpatterns+=static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
    urlpatterns+=static(settings.MEDIA_URL,document_root=settings.MEDIA_ROOT)