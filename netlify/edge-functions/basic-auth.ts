import { Context } from "https://edge.netlify.com";

export default async (request: Request, context: Context) => {
    const USERNAME = Deno.env.get("BASIC_USERNAME");
    const PASSWORD = Deno.env.get("BASIC_PASSWORD");
    
    const auth = request.headers.get("authorization");
    if (!auth || !auth.startsWith("Basic ")) return unauthorized();
    
    const [user, pass] = atob(auth.split(" ")[1]).split(":");
    
    if (user !== USERNAME || pass !== PASSWORD) return unauthorized();
    
    return context.next();
};

function unauthorized() {
    return new Response("🔐 Unauthorized", {
        status: 401,
        headers: {
            "WWW-Authenticate": 'Basic realm="Protected"',
            "Content-Type": "text/plain"
        }
    });
}
