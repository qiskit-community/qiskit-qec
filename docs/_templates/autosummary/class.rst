{#
   The general principle of this is that we manually document methods and attributes here in
   the same file.  By default, we document
   all methods, including those defined by parent classes.
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
{#
   The documentation from `__init__` is automatically included in the
   standard class documentation, so we don't want to repeat it.
-#}
{% block methods_summary %}{% set wanted_methods = (methods | reject('==', '__init__') | list) %}{% if wanted_methods %}
   .. rubric:: Methods

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
