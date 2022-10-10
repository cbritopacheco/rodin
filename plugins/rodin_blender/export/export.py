import bpy

class RODIN_OT_Export(bpy.types.Operator):
    bl_idname = 'rodin.export'
    bl_label = 'RodinExport'
    bl_description = 'Exports the mesh into a Rodin supported format.'

    def __init__(self):
        self. writepath = None

    def execute(self, context):
        scene = context.scene
        for obj in scene.objects:
            if obj is not None and obj.type == 'MESH':
                mesh = obj



