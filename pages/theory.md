---
layout: page
show_meta: false
title: "Theories"
subheadline: "Theories"
header:
   image_fullwidth: "head_beckman.jpg"
image:
    title: theory_title.png
permalink: "/theory/"
---

Here are all the categories of learning quantum chemistry, condensed matter physics and quantum computation.

<ul>
    {% for post in site.categories.theory %}
    <li><a href="{{ site.url }}{{ site.baseurl }}{{ post.url }}">{{ post.title }}</a></li>
    {% endfor %}
</ul>