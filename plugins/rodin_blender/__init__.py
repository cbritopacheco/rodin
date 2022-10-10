import bpy

import rodin_blender.gui as gui

bl_info = {
    'name' : 'Rodin',
    'author' : 'Carlos Brito-Pacheco',
    'description' : 'Rodin\'s companion tool for use with Blender.',
    'blender' : (3, 0, 0),
    'version' : (0, 0, 1),
    'location' : 'View3D',
    'warning' : '',
    'wiki_url' : 'https://cbritopacheco.github.io/rodin/',
    'tracker_url' : 'https://github.com/cbritopacheco/rodin',
    'category' : 'Mesh'
}

def register():
    gui.register()

def unregister():
    gui.unregister()

if __name__ == "__main__":
    register()

