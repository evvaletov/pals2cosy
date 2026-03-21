project = 'pals2cosy'
copyright = '2026, Eremey Valetov'
author = 'Eremey Valetov'
release = '0.4.0'

extensions = ['myst_parser']
myst_enable_extensions = ['dollarmath', 'amsmath']

source_suffix = {'.md': 'markdown'}
master_doc = 'index'
exclude_patterns = ['_build']

html_theme = 'alabaster'
html_theme_options = {
    'description': 'PALS lattice → COSY INFINITY converter',
    'fixed_sidebar': True,
}

latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '11pt',
}
latex_documents = [
    (master_doc, 'pals2cosy.tex', 'pals2cosy Documentation',
     'Eremey Valetov', 'manual'),
]
