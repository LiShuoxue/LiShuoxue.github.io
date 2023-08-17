---
layout: page
show_meta: false
title: "Theories"
subheadline: "Theories"
header:
   image_fullwidth: "head_goldgate.jpg"
permalink: "/theory/"
---
<ul>
    {% for post in site.categories.theory %}
    <li><a href="{{ site.url }}{{ site.baseurl }}{{ post.url }}">{{ post.title }}</a></li>
    {% endfor %}
</ul>