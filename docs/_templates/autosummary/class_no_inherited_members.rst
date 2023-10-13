{#
   This is very similar to the default class template, except this one is used
   when we don't want to generate any inherited methods.
-#}

{{ objname | escape | underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
{#-
   Avoid having autodoc populate the class with the members we're about to
   summarize to avoid duplication.
#}
   :no-members:
   :show-inheritance:
   :no-inherited-members:
   :no-special-members:

{% block methods_summary %}{% set wanted_methods = (methods | reject('in', inherited_members) | reject('==', '__init__') | list) %}{% if wanted_methods %}
   .. rubric:: Methods Defined Here

{% for item in wanted_methods %}
   .. automethod:: {{ name }}.{{ item }}
{%- endfor %}
{% endif %}{% endblock %}

{% block attributes_summary %}{% if attributes %}
   .. rubric:: Attributes
{# Attributes should all be summarized directly on the same page. -#}
{% for item in attributes %}
   .. autoattribute:: {{ item }}
{%- endfor %}
{% endif %}{% endblock -%}
